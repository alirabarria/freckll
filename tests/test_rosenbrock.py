import logging

import numpy as np

from freckll.kinetics import AltitudeSolveError
from freckll.solver.rosenbrock import Rosenbrock, update_timestep
from freckll.solver.transform import LogTransform, UnityTransform


def test_altitude_error_during_candidate_validation_rejects_step(monkeypatch, caplog):
    """A failed candidate validation should reduce the step and retry."""
    attempted_steps = []

    def fake_step(f, jac, y, t, h):
        attempted_steps.append(h)
        return y + 0.1, np.full_like(y, 1e-4)

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

    with caplog.at_level(logging.INFO):
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
            trace_attempts=True,
        )

    assert result["success"] is True
    assert validation_calls == 2
    assert attempted_steps == [0.5, 0.05]
    assert result["times"][-1] == 0.05
    messages = [record.getMessage() for record in caplog.records]
    assert any(
        "Rosenbrock attempt:" in message
        and "h=5.00000000000000000E-01" in message
        for message in messages
    )
    assert any(
        "reason=altitude_during_candidate_validation" in message
        for message in messages
    )
    assert any(
        "Rosenbrock accepted:" in message
        and "h=5.00000000000000028E-02" in message
        for message in messages
    )


def test_failed_log_solve_evaluates_diagnostics_in_transformed_space(monkeypatch):
    """Failure diagnostics must not apply LogTransform.inverse twice."""
    physical_y = np.array([2.0])
    transformed_y = np.log(physical_y)
    f_inputs = []
    jac_inputs = []

    def fake_step(f, jac, y, t, h):
        return y.copy(), np.full_like(y, 1e-4)

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


def test_update_timestep_bounds_growth_and_reduction():
    """The controller must remain finite and bound each multiplicative change."""
    assert update_timestep(1.0, 1e-3, 0.0) == 2.0
    assert update_timestep(1.0, 1e-3, 1e3) == 0.5


def test_local_error_rejects_candidate_before_advancing(monkeypatch, caplog):
    """A candidate with delta above rtol must be retried at the same time."""
    attempted_steps = []
    attempted_times = []
    errors = iter((1e-2, 1e-4))

    def fake_step(f, jac, y, t, h):
        attempted_steps.append(h)
        attempted_times.append(t)
        return y + 0.1, np.full_like(y, next(errors))

    monkeypatch.setattr(
        "freckll.solver.rosenbrock.step_second_order_rosenbrock",
        fake_step,
    )

    solver = Rosenbrock.__new__(Rosenbrock)
    solver._logger = logging.getLogger("freckll.test_rosenbrock")

    with caplog.at_level(logging.INFO):
        result = solver._run_solver(
            f=lambda t, y: np.zeros_like(y),
            jac=lambda t, y: np.eye(y.size),
            y0=np.array([1.0]),
            t0=0.0,
            t1=0.05,
            num_species=1,
            transform=UnityTransform(),
            initial_step=0.1,
            rtol=1e-3,
            maxiter=3,
            trace_attempts=True,
        )

    assert result["success"] is True
    assert attempted_steps == [0.1, 0.05]
    assert attempted_times == [0.0, 0.0]
    assert result["times"][-1] == 0.05
    assert any(
        "reason=local_error" in record.getMessage()
        and "delta=1.00000000000000002E-02" in record.getMessage()
        for record in caplog.records
    )
