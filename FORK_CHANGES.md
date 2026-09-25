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

## 5. Restore local-error rejection and bound Rosenbrock timestep changes

**Commit:** `Restore Rosenbrock local-error rejection` (development change on
`codex/rosenbrock-step-recovery`).

### Problem

FRECKLL's Rosenbrock step returns a second-order candidate together with an
embedded error estimate: the absolute difference between its second-order and
first-order solutions. The solver reduced that array to the scalar `delta`,
but did so only after it had already accepted the candidate, copied it into
the integration state, and advanced the simulation time. Consequently,
`delta > rtol` could only reduce the *next* timestep; it could not reject the
inaccurate candidate that produced the large error.

This was visible in a matched W39/Veillet-2024 replay on Geryon. One finite
candidate had approximately `delta = 26.47` for `rtol = 1e-3` and a maximum
abundance of `1.19e5`. FRECKLL accepted that state and reduced only the next
timestep. Once the accepted state had been corrupted, later error estimates
could become exactly zero and the original controller divided by zero,
allowing an infinite proposed timestep.

The behavior also differed from VULCAN, on which this Rosenbrock
implementation is based. VULCAN accepts a step only when its local error is
within tolerance, rejects and retries otherwise, substitutes a small finite
error when `delta == 0`, and bounds the multiplicative timestep change.

### Change

FRECKLL now performs the operations in this order:

1. build and physically validate the candidate state;
2. apply the existing `atol` abundance mask to the embedded error;
3. compute the existing scalar `delta`;
4. calculate a bounded candidate for the next timestep;
5. reject the candidate without changing `y` or `t` when `delta > rtol`; and
6. accept and advance only when the candidate passes the local-error test.

The default timestep-factor bounds are `0.5` and `2.0`. These are per-step
ratio limits, not absolute timestep limits: repeated accepted steps can still
span the many orders of magnitude required by atmospheric chemistry. An
exactly zero `delta` is treated as `0.01 * rtol` for timestep selection, which
results in the bounded maximum growth rather than infinity.

With attempt tracing enabled, these retries are logged as
`reason=local_error`, including `t`, `h`, `delta`, `rtol`, the proposed next
`h`, and the candidate range.

### What this change deliberately does not alter

This is a control-flow correction, not a redefinition of FRECKLL's error
model. In particular:

- `atol` still acts as the abundance threshold that excludes species from the
  scalar error comparison;
- under `LogTransform`, the embedded error remains a difference in
  log-abundance coordinates rather than a conventional
  `ATOL + RTOL * abs(y)` weighted norm;
- the convergence criteria `df_criteria` and `dfdt_criteria` are unchanged;
  and
- this cannot repair an inaccurate derivative or Jacobian caused by
  catastrophic cancellation. It is intended to prevent an obviously poor
  embedded Rosenbrock candidate from being accepted after that numerical
  problem has manifested.

### Tests

`tests/test_rosenbrock.py` now verifies that:

- a candidate with `delta > rtol` is retried at the same simulation time;
- the rejected candidate never becomes the accepted state;
- the retry uses the bounded reduction factor;
- `delta == 0` produces finite growth capped at a factor of two; and
- a very large error produces a reduction capped at a factor of one half.

### Scientific and campaign implications

This change can alter every Rosenbrock trajectory whose historical run
accepted at least one candidate with `delta > rtol`, so patched campaigns
must use a new output directory and record the FRECKLL commit. Comparisons
against old campaigns are diagnostic comparisons, not bitwise continuations.
The W39 Monte Carlo campaign should be restarted only after a small matched
seed replay confirms that the new local-error rejections occur as expected
and that the solver still reaches the intended physical solution.

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
