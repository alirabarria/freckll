# Changes maintained in this fork

This document records changes made in Arturo Lira-Barria's FRECKLL fork for
the FETUKINES Monte Carlo analysis. Its purpose is to make every numerical or
reproducibility change explicit, testable, and citable. Each implementation
change is kept in a separate Git commit whenever practical.

The current development branch is `codex/rosenbrock-step-recovery`. These
changes have not yet been merged into this fork's `master` branch.

## 1. Recover from altitude failures during Rosenbrock candidate validation

**Commit:** `1e3dff7` (`Recover from altitude errors during Rosenbrock validation`)

### Problem

The Rosenbrock solver already caught `AltitudeSolveError` while constructing a
candidate step. It then evaluated the derivative once more at the proposed
state before accepting the step:

```python
test_f = f(t, transform.transform(y_new))
```

That validation call was outside the existing exception handler. A candidate
state could therefore trigger `AltitudeSolveError` during validation and abort
the complete atmospheric realization, even though rejecting the candidate and
retrying with a smaller timestep was possible.

This failure was observed in W39/Veillet-2024 Monte Carlo runs.

### Change

The validation call now catches `AltitudeSolveError`. On failure, FRECKLL:

1. leaves the last accepted state unchanged;
2. rejects the proposed state;
3. multiplies the timestep by `timestep_reject_factor`;
4. stops if the reduced timestep is below `minimum_step`; and
5. otherwise retries the Rosenbrock step.

This is a defensive recovery change. It does not alter a step that passes
candidate validation.

### Test

`tests/test_rosenbrock.py` constructs a deterministic candidate that raises
`AltitudeSolveError` on its first validation. The test verifies that the first
step is rejected, the timestep is reduced, and the subsequent retry succeeds.

### Remaining limitations

This commit does not yet address non-finite candidate states, overflow in
`LogTransform`, failure diagnostics, integration-range overshoot, or the
meaning of the solver's `success` flag. Those concerns must be handled and
tested separately.

## 2. Load Venot reaction files in deterministic order

**Commit:** `15ff333` (`Load Venot reactions in deterministic order`)

### Problem

The Venot network loader iterated directly over:

```python
directory.glob("*.dat")
```

Filesystem iteration order is not guaranteed. For the same Veillet-2024
network, macOS/APFS and the Geryon Linux filesystem returned the same reaction
files in different orders. FRECKLL consequently assembled the same 1,310
reactions in different sequences.

FETUKINES generates a deterministic series of Monte Carlo draws from each
seed and assigns those draws to `network.reaction_calls` by position. Before
this fix, the seed reproduced the random numbers but not the reaction to which
each number was assigned. A direct comparison found:

- identical multisets of all 1,310 reactions and coefficients;
- identical initial random-number sequences;
- different row order for 1,245 of the 1,310 random-draw assignments; and
- different final perturbation factors for 989 reactions.

Thus, a seed and run ID were not sufficient to reproduce a Monte Carlo
realization across filesystems.

### Change

Reaction files are now sorted before they are parsed:

```python
reaction_files = sorted(directory.glob("*.dat"))
```

The species list produced by `infer_composition` is also sorted by its string
representation instead of being returned directly from a Python `set`. This
removes a second source of platform- and hash-dependent ordering.

### Test

`tests/test_venot.py::test_load_reactions_sorts_input_files` deliberately
supplies reaction files in reverse order and verifies that `load_reactions`
returns their reaction calls in canonical filename order.

### Scientific and campaign implications

The chemical network itself is unchanged; only its in-memory order becomes
deterministic. Nevertheless, this changes the mapping between historical
`run_id`/seed pairs and reactions for campaigns produced with the unsorted
loader. A campaign started with an older FRECKLL version must not be continued
in the same output directory after adopting this commit.

New campaigns should record at least:

- the FRECKLL Git commit;
- the FETUKINES Git commit;
- the explicit base seed;
- the reaction-network identity or checksums; and
- the generated per-run perturbation table.

## 3. Evaluate unsuccessful-solve diagnostics in transformed space

**Commit:** `aafb41f` (`Evaluate failed-solve diagnostics in transformed space`)

### Problem

The Rosenbrock integration loop keeps its accepted state `y` in physical
abundance space. A candidate is produced in transformed coordinates and then
returned to physical coordinates before it is accepted:

```python
y_new = transform.inverse_transform(y_new)
y = np.copy(y_new)
```

The derivative and Jacobian callables, however, expect their input in the
selected transformed space because both invert that transform internally.
Most calls in the integration loop correctly use a transformed state, for
example:

```python
test_f = f(t, transform.transform(y_new))
```

When an integration ended unsuccessfully, FRECKLL assembled additional
diagnostics with the physical state directly:

```python
extra["dndt"] = f(t, y)
extra["jac"] = jac(t, y)
```

For `LogTransform`, this made `f` and `jac` interpret physical abundances as
log-abundances and apply `exp(y)` to them. This is not a pending transform from
the preceding step: `y` is already the physical state. The mismatch could
produce overflow, non-finite chemistry, or an `AltitudeSolveError` while
reporting an otherwise ordinary unsuccessful termination such as `maxiter`.
The secondary exception could therefore hide the actual reason the solver
stopped.

### Change

Failure diagnostics now receive the same transformed representation used by
the rest of the solver:

```python
diagnostic_y = transform.transform(y)
extra["dndt"] = f(t, diagnostic_y)
extra["jac"] = jac(t, diagnostic_y)
```

This change does not modify accepted steps, timestep selection, convergence
criteria, or successful integrations. It only corrects diagnostic evaluation
after `success=False`.

### Test

`tests/test_rosenbrock.py::test_failed_log_solve_evaluates_diagnostics_in_transformed_space`
forces a deterministic `LogTransform` integration to terminate at `maxiter`.
It verifies that both final diagnostic callables receive `log(y)`, rather than
the already physical `y`.

### Scientific and campaign implications

Successful historical runs are numerically unaffected. Some historical runs
reported as worker exceptions may instead have completed the integration loop
with `success=False` and then crashed only while assembling failure
diagnostics. Those runs should not be reclassified as converged, but the patch
allows their genuine termination state and diagnostic arrays to be retained.

## 4. Optional trace of every Rosenbrock timestep attempt

**Commit:** `Trace every Rosenbrock timestep attempt` (temporary diagnostic
change on `codex/rosenbrock-step-recovery`).

### Problem

FRECKLL normally logs accepted steps but does not identify every attempted
timestep or the exact branch that rejected a candidate. Consequently, two
platforms can show different first accepted timesteps without revealing where
their numerical trajectories first diverged.

### Change

The Rosenbrock solver accepts an opt-in `trace_attempts` argument. The same
mode can be enabled without changing a calling application by setting
`FRECKLL_TRACE_ROSENBROCK_ATTEMPTS=1`. When enabled, it logs each attempted
`iteration`, `t`, and `h`. It then records whether the candidate was accepted
or rejected. Rejections distinguish:

- an `AltitudeSolveError` during the Rosenbrock stages;
- a stage derivative containing `NaN`;
- an `AltitudeSolveError` during candidate validation; and
- a candidate with non-finite/negative abundances or a `NaN` validation
  derivative.

Accepted-step records include the estimated absolute transformed-coordinate
error `delta`, the next timestep, candidate range, error range, validation
derivative range, and count of infinite validation derivatives.

Tracing is disabled by default and does not change the solver decision logic.
It is intended only for matched cross-platform replays because it can produce
large log files.

### Test

The candidate-altitude regression test now enables tracing and verifies that
the rejected initial attempt and accepted retry are both recorded with their
respective timesteps and rejection reason.

## Reproducibility policy for future changes

Future modifications should be added here with:

1. the commit identifier;
2. the observed failure or ambiguity;
3. the exact behavioral change;
4. the associated regression test;
5. whether existing campaigns remain comparable; and
6. any required migration or rerun.

Numerical experiments that deliberately change the integration method, such
as replacing FRECKLL's current Rosenbrock error controller with a conventional
`ATOL + RTOL * abs(y)` scale, should be developed on a separate experimental
branch and must not be described as compatibility fixes.
