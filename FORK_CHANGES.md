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
