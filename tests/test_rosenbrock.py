import logging

import numpy as np

from freckll.kinetics import AltitudeSolveError
from freckll.solver.rosenbrock import Rosenbrock
from freckll.solver.transform import LogTransform, UnityTransform


def test_altitude_error_during_candidate_validation_rejects_step(monkeypatch):
    """A failed candidate validation should reduce the step and retry."""
    attempted_steps = []

    def fake_step(f, jac, y, t, h):
        attempted_steps.append(h)
        return y + 0.1, np.full_like(y, 0.01)

    monkeypatch.setattr(
        "freckll.solver.rosenbrock.step_second_order_rosenbrock",
        fake_step,
    )

    validation_calls = 0

    def f(t, y):
        nonlocal validation_calls
        validation_calls += 1
        if validation_calls == 1:
            raise AltitudeSolveError
        return np.zeros_like(y)

    solver = Rosenbrock.__new__(Rosenbrock)
    solver._logger = logging.getLogger("freckll.test_rosenbrock")

    result = solver._run_solver(
        f=f,
        jac=lambda t, y: np.eye(y.size),
        y0=np.array([1.0]),
        t0=0.0,
        t1=0.05,
        num_species=1,
        transform=UnityTransform(),
        initial_step=0.5,
        timestep_reject_factor=0.1,
        minimum_step=1e-6,
        maxiter=5,
    )

    assert result["success"] is True
    assert validation_calls == 2
    assert attempted_steps == [0.5, 0.05]
    assert result["times"][-1] == 0.05


def test_failed_log_solve_evaluates_diagnostics_in_transformed_space(monkeypatch):
    """Failure diagnostics must not apply LogTransform.inverse twice."""
    physical_y = np.array([2.0])
    transformed_y = np.log(physical_y)
    f_inputs = []
    jac_inputs = []

    def fake_step(f, jac, y, t, h):
        return y.copy(), np.full_like(y, 0.01)

    monkeypatch.setattr(
        "freckll.solver.rosenbrock.step_second_order_rosenbrock",
        fake_step,
    )

    def f(t, y):
        f_inputs.append(y.copy())
        return np.zeros_like(y)

    def jac(t, y):
        jac_inputs.append(y.copy())
        return np.eye(y.size)

    solver = Rosenbrock.__new__(Rosenbrock)
    solver._logger = logging.getLogger("freckll.test_rosenbrock")
    result = solver._run_solver(
        f=f,
        jac=jac,
        y0=physical_y,
        t0=0.0,
        t1=1.0,
        num_species=1,
        transform=LogTransform(),
        initial_step=0.1,
        maxiter=1,
        df_criteria=-1.0,
        dfdt_criteria=-1.0,
    )

    assert result["success"] is False
    assert len(f_inputs) == 2  # candidate validation and failure diagnostic
    np.testing.assert_allclose(f_inputs[-1], transformed_y)
    assert len(jac_inputs) == 1
    np.testing.assert_allclose(jac_inputs[-1], transformed_y)
