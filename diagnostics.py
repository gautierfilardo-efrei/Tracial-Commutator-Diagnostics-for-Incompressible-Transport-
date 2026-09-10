"""Transport diagnostics on the unit torus. No Navier--Stokes blow-up detector."""
from itertools import product
from operator import index
import numpy as np


def _positive_integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a positive integer")
    try:
        value = index(value)
    except TypeError as error:
        raise ValueError(f"{name} must be a positive integer") from error
    if value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _points(values, dimension=None):
    points = np.asarray(values, dtype=float)
    if points.ndim != 2 or points.shape[1] < 1:
        raise ValueError("points must have shape (samples, positive dimension)")
    if dimension is not None and points.shape[1] != dimension:
        raise ValueError(f"points must have exactly {dimension} coordinates")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must be finite")
    return points


def mode_count(k, dimension=3):
    k = _positive_integer(k, "k")
    dimension = _positive_integer(dimension, "dimension")
    return (2 * k + 1) ** dimension - 1


def cube_modes(k, dimension=3):
    k = _positive_integer(k, "k")
    dimension = _positive_integer(dimension, "dimension")
    return np.array([v for v in product(range(-k, k + 1), repeat=dimension)
                     if any(v)], dtype=float)


def energy_coefficient(k, dimension=3):
    """Sum of k_1^2 over the nonzero Fourier cube."""
    k = _positive_integer(k, "k")
    dimension = _positive_integer(dimension, "dimension")
    return k * (k + 1) * (2 * k + 1) ** dimension // 3


def universal_bound(k, dimension=3):
    return float(np.log1p(4 * mode_count(k, dimension)))


def cube_increment_direct(displacement, k, block_size=128):
    """Stable reference sum of 4 sin^2(pi k.delta), O(N K^d)."""
    delta = _points(displacement)
    delta = delta - np.rint(delta)
    block_size = _positive_integer(block_size, "block_size")
    modes = cube_modes(k, delta.shape[1])
    result = np.zeros(len(delta))
    for start in range(0, len(modes), block_size):
        phase = delta @ modes[start:start + block_size].T
        result += 4 * np.sum(np.sin(np.pi * phase) ** 2, axis=1)
    return result


def cube_increment(displacement, k):
    """Dirichlet-kernel formula, with cancellation-safe direct fallback."""
    delta = _points(displacement)
    k = _positive_integer(k, "k")
    mode_count(k, delta.shape[1])
    delta = delta - np.rint(delta)
    width = 2 * k + 1
    full_count = width ** delta.shape[1]
    kernels = width * np.sinc(width * delta) / np.sinc(delta)
    result = 2 * (full_count - np.prod(kernels, axis=1))
    sensitive = result < 1e-7 * full_count
    if np.any(sensitive):
        result[sensitive] = cube_increment_direct(delta[sensitive], k)
    upper = 4 * mode_count(k, delta.shape[1])
    if np.any(result < -1e-9) or np.any(result > upper + 1e-7 * full_count):
        raise ArithmeticError("Fourier increment violates its analytical bound")
    return result


def cube_integrand(displacement, k, normalized=False):
    delta = np.asarray(displacement)
    h = cube_increment(delta, k)
    if normalized:
        h = h / mode_count(k, delta.shape[1])
    return np.log1p(h)


def mean_se(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.all(np.isfinite(values)):
        raise ValueError("a finite one-dimensional array of at least two samples is required")
    return float(values.mean()), float(values.std(ddof=1) / np.sqrt(len(values)))


def abc_velocity(x):
    """A=B=C=1, period one; curl(v)=2*pi*v and Laplacian(v)=-4*pi^2*v."""
    z = 2 * np.pi * _points(x, dimension=3)
    return np.column_stack((np.sin(z[:, 2]) + np.cos(z[:, 1]),
                            np.sin(z[:, 0]) + np.cos(z[:, 2]),
                            np.sin(z[:, 1]) + np.cos(z[:, 0])))


def abc_flow(x, lag, steps=128, viscosity=0.01, start_time=0.0):
    """Lifted trajectories of u(x,t)=exp(-4*pi^2*nu*t)*v(x), via RK4.

    This explicit velocity solves unforced periodic Navier--Stokes with
    p=-exp(-8*pi^2*nu*t)*|v|^2/2. Only the particle ODE is discretized.
    """
    steps = _positive_integer(steps, "steps")
    if not np.all(np.isfinite([lag, viscosity, start_time])) or lag < 0 or viscosity < 0:
        raise ValueError("invalid flow parameters")
    c = 4 * np.pi ** 2 * viscosity
    clock = lag if c == 0 else np.exp(-c * start_time) * (-np.expm1(-c * lag)) / c
    y = _points(x, dimension=3).copy()
    dt = clock / steps
    for _ in range(steps):
        k1 = abc_velocity(y)
        k2 = abc_velocity(y + dt * k1 / 2)
        k3 = abc_velocity(y + dt * k2 / 2)
        k4 = abc_velocity(y + dt * k3)
        y += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return y
