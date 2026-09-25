# FRECKLL

[![Release](https://img.shields.io/github/v/release/ahmed-f-alrefaie/freckll)](https://img.shields.io/github/v/release/ahmed-f-alrefaie/freckll)
[![Build status](https://img.shields.io/github/actions/workflow/status/ahmed-f-alrefaie/freckll/main.yml?branch=main)](https://github.com/ahmed-f-alrefaie/freckll/actions/workflows/main.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/ahmed-f-alrefaie/freckll/branch/main/graph/badge.svg)](https://codecov.io/gh/ahmed-f-alrefaie/freckll)
[![Commit activity](https://img.shields.io/github/commit-activity/m/ahmed-f-alrefaie/freckll)](https://img.shields.io/github/commit-activity/m/ahmed-f-alrefaie/freckll)
[![License](https://img.shields.io/github/license/ahmed-f-alrefaie/freckll)](https://img.shields.io/github/license/ahmed-f-alrefaie/freckll)

Fast Disequilibrium chemistry

- **Github repository**: <https://github.com/ahmed-f-alrefaie/freckll/>
- **Documentation** <https://ahmed-f-alrefaie.github.io/freckll/>
- **Changes maintained in this fork**: [FORK_CHANGES.md](FORK_CHANGES.md).
  This changelog documents each numerical and reproducibility fix by commit,
  including regression tests and implications for existing campaigns and pull
  requests.

This fork restores rejection of Rosenbrock candidates whose embedded local
error exceeds `rtol`, bounds each multiplicative timestep change, and prevents
an exactly zero error estimate from producing an infinite timestep. The
detailed changelog explains the numerical scope and limitations of this
control-flow correction; it does not redefine FRECKLL's `atol`/`rtol` metric
or claim to fix inaccurate chemistry derivatives and Jacobians.

## Test case

```bash
cd freckll
pip install .
```

```bash
cd freckll/examples/inputs

freckll hd209458_full.yml -o result_full.h5 --overwrite --plot
```
