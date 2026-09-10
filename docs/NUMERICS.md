# Numerical conventions and limitations

- Unit torus: all coordinates have period one. Fourier modes are
  `exp(2*pi*i*k.x)` and all logarithms are natural.
- The Fourier list is the complete nonzero cube `||k||_infinity <= K`.
  It is not a cutoff applied to the velocity field.
- Monte Carlo sample: `N=4096`, `numpy.random.default_rng(20260910)`.
  The sample coordinates are provided as `results/sample_points.npy`.
- Standard errors use sample standard deviation (`ddof=1`) divided by
  `sqrt(N)`. They measure sampling uncertainty only.
- The translation and cat-map reference values follow exact integral formulas.
- The ABC field is `u(x,t)=exp(-4*pi*pi*nu*t)*v(x)`, with A=B=C=1 and
  `nu=0.01`. Its pressure is `-exp(-8*pi*pi*nu*t)*|v|**2/2`.
- RK4 integrates the autonomous field `v` over the exact effective time
  `exp(-c*t)*(-expm1(-c*h))/c`, with `c=4*pi*pi*nu` and the continuous
  limit `h` at `nu=0`. Lifted positions are returned.
- Small-time checks use `K=4`, 16 steps per lag and lags from `2e-3` to
  `2e-5`. The exact spatial energy limit and the same-sample limit are distinct.
- The lag-one refinement experiment compares 64, 128 and 256 RK4 steps.
  This is an empirical check, not a certified trajectory-error bound.
- Shears use `h=0.2`, `K=4`, a deterministic grid of `2**15` midpoints and
  frequencies `1,3,9,27,81`. No Monte Carlo SE is assigned to this quadrature.
- The Dirichlet product uses a direct nonnegative sine-square fallback when
  its increment is less than `1e-7*(2*K+1)**d`, to control cancellation.
- Frequency cutoffs, dimensions and step counts must be positive integers.
  Nonfinite input arrays are rejected rather than silently propagated.

The original table audit transcribes nine supplied values and compares them
with the universal bound. It is a record of contradictions in that table,
not a newly simulated fluid dataset.
