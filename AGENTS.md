# AGENTS.md

## Project overview

This repository is for C++ numerical experiments on a stochastic path-following model near a simple parameterized curve.

The mathematical draft is about a particle/person moving near a target curve `Gamma(t)` while trying to both return to the curve and advance in the positive parameter direction. The simulation should generate true motion points and optionally sparse, noisy observed points.

The core model uses:

- A target parameter curve `Gamma(t)`.
- A locally tracked nearest parameter `t_k` and nearest point `Q_k = Gamma(t_k)`.
- A forward lookahead point `Q_+ = Gamma(t_k + T_k)` determined by arc length `L`.
- An expected direction

  ```text
  r_exp(P_k, L) = normalize(Q_+ - P_k)
  ```

- Direction dynamics with turning inertia and angular noise.
- Speed dynamics with speed inertia and speed noise.
- Optional sparse sampling and measurement noise.

The code should be written as a simulation/experiment framework, not just a one-off script.

## Important mathematical conventions

Use radians everywhere.

Use the standard 2D position convention:

```text
P = (x, y)
r(theta) = (cos(theta), sin(theta))
```

For the standard ellipse, use:

```text
Gamma(t) = (a cos t, b sin t)
Gamma'(t) = (-a sin t, b cos t)
|Gamma'(t)| = sqrt(a^2 sin^2 t + b^2 cos^2 t)
```

For a straight vertical line example, use:

```text
Gamma(t) = (0, t)
positive direction = +y
Q_+ = (0, y + L)
r_exp(P, L) = normalize((-x, L))
theta_exp = atan2(L, -x)
```

Do not replace `atan2(y, x)` with `atan(y/x)`. The quadrant matters.

For closed curves, prefer an unfolded real-valued parameter `t in R` during simulation and only apply periodic wrapping when evaluating `Gamma(t)`. This avoids jumps at `0` and `2*pi`.

## Core simulation state

Maintain a state similar to:

```cpp
struct State {
    Vec2 P;        // true position
    double theta;  // current heading angle, radians
    double v;      // current speed
    double t;      // locally tracked nearest curve parameter
};
```

A full simulation step should roughly follow this order:

1. Given the current `P_k`, `theta_k`, `v_k`, and tracked parameter `t_k`, compute the expected direction.
2. Use `t_k` to compute `Q_k` and the lookahead point `Q_{+,k}`.
3. Compute `theta_exp_{k+1} = atan2(r_exp.y, r_exp.x)`.
4. Update heading:

   ```text
   theta_{k+1} = theta_k
               + alpha * wrap(theta_exp_{k+1} - theta_k)
               + epsilon_theta_{k+1}
   ```

5. Update speed:

   ```text
   v_{k+1} = v_k + beta * (v_exp - v_k) + epsilon_v_{k+1}
   ```

   If the configuration requests nonnegative speeds, clamp with:

   ```text
   v_{k+1} = max(0, v_{k+1})
   ```

6. Update true position:

   ```text
   P_{k+1} = P_k + v_{k+1} * dt * (cos(theta_{k+1}), sin(theta_{k+1}))
   ```

7. Update the locally tracked nearest parameter `t_{k+1}` by prediction plus local projection.
8. If `k + 1` is a sampling step, write the true point and an optional measured point.

## Noise model

Use reproducible random generation.

Recommended default implementation:

```cpp
std::mt19937_64 rng(seed);
std::normal_distribution<double> angle_noise(0.0, sigma_theta);
std::normal_distribution<double> speed_noise(0.0, sigma_v);
std::normal_distribution<double> measure_noise(0.0, sigma_measure);
```

Noise terms:

- `epsilon_theta`: angular control noise.
- `epsilon_v`: speed control noise.
- `(delta_x, delta_y)`: measurement noise for observed positions.

The measurement model should be:

```text
P_scan = P_true + (delta_x, delta_y)
delta_x, delta_y iid N(0, sigma_measure^2)
```

Treat movement noise and measurement noise as separate concepts. Do not apply measurement noise to the true state update.

Optional later extension: a von Mises distribution may replace Gaussian angular noise for large circular angular uncertainty. Do not implement it unless explicitly requested.

## Local nearest-parameter tracking

Do not recompute the nearest point by global search at every step unless explicitly requested. This is slow and can incorrectly jump to a different nearby branch of the curve.

Track `t_k` as a state variable.

Prediction:

```text
t_{k+1}^{(0)} = t_k + ((P_{k+1} - P_k) dot Gamma'(t_k)) / |Gamma'(t_k)|^2
```

Then define a local interval:

```text
I_k = [t_{k+1}^{(0)} - delta, t_{k+1}^{(0)} + delta]
```

Within this interval, minimize:

```text
h(t) = |P_{k+1} - Gamma(t)|^2
```

A Newton candidate may be found by solving the stationary condition:

```text
(P_{k+1} - Gamma(t)) dot Gamma'(t) = 0
```

Newton iteration:

```text
t_next = t - ((P - Gamma(t)) dot Gamma'(t)) /
             (-|Gamma'(t)|^2 + (P - Gamma(t)) dot Gamma''(t))
```

Implementation requirements:

- Limit Newton iterations with a configurable `max_newton_iters`.
- Stop when `abs(t_next - t) < newton_tol`.
- If the Newton point converges and lies in `I_k`, compare it with both interval endpoints.
- If Newton fails or leaves `I_k`, choose the better endpoint.
- The final `t_{k+1}` must be the candidate with the smallest `h(t)` among the valid candidates.

This local projection rule is important for preventing branch jumps on curves with spatially nearby but parametrically distant segments.

## Lookahead arc-length computation

The lookahead parameter increment `T` is defined by:

```text
integral from t_k to t_k + T of |Gamma'(t)| dt = L
```

Provide multiple strategies when practical:

1. Exact/simple strategy for simple curves.
2. First-order Taylor strategy:

   ```text
   T ≈ L / |Gamma'(t_k)|
   ```

3. Higher-order Taylor/Lagrange-inversion strategy when derivatives are available.
4. Robust numerical fallback, e.g. bisection on accumulated arc length or adaptive Simpson, if requested.

For the standard ellipse, at least implement first-order and second-order approximations consistently with:

```text
f(t) = sqrt(b^2 + (a^2 - b^2) sin^2 t)
```

and

```text
T ≈ L / f(t_0) - ((a^2 - b^2) sin(2 t_0)) / (4 f(t_0)^4) * L^2
```

Be careful: the denominator of the second-order term is `f(t_0)^4`, not `f(t_0)^2`.

## Expected code organization

Prefer a small, readable C++ project structure such as:

```text
include/
  vec2.hpp
  curve.hpp
  simulation.hpp
  rng.hpp
src/
  main.cpp
  simulation.cpp
  curves.cpp
tests/
  test_geometry.cpp
  test_line_trace.cpp
  test_ellipse_trace.cpp
outputs/
  *.csv
CMakeLists.txt
AGENTS.md
README.md
```

If the repository already has a different structure, adapt to it instead of forcing this layout.

## C++ style requirements

- Use modern C++17 or C++20.
- Prefer `double` for numerical state and geometry.
- Keep geometry utilities small and explicit.
- Avoid hidden global mutable state except a deliberate RNG object owned by the simulation/config.
- Make randomness reproducible through a seed in the configuration.
- Prefer clear names over terse names.
- Keep formulas close to the mathematical notation used in the draft.
- Use `std::atan2`, not plain `atan`, for angles.
- Implement `wrap_angle` to map angles to `(-pi, pi]` or a clearly documented equivalent interval.
- Avoid unnecessary external dependencies. Use the standard library unless a dependency is clearly justified.

## Configuration

Prefer a configurable experiment interface. It can be command-line flags, a small config struct, or a simple config file.

Important parameters:

```text
curve_type
n_steps
dt
lookahead_L
alpha
beta
v_exp
sigma_theta
sigma_v
sigma_measure
sample_interval_N
seed
local_projection_delta
max_newton_iters
newton_tol
speed_nonnegative
arc_length_method
output_csv_path
```

## CSV output requirements

Write enough columns to debug and reproduce experiments.

Recommended columns:

```text
index,
time,
true_x,
true_y,
scan_x,
scan_y,
theta,
theta_exp,
v,
t_min,
T,
epsilon_theta,
epsilon_v,
measure_dx,
measure_dy,
was_sampled
```

For comparison experiments, include enough metadata at the top of the file or in a sidecar file to record parameters. If metadata lines are placed in the CSV, make sure the plotting/reading utilities can skip them reliably.

## Validation and tests

Add small deterministic tests wherever possible.

Minimum useful tests:

1. `wrap_angle` maps results into the documented interval.
2. `atan2` direction for the line example: for `P=(0,y)`, expected heading is `pi/2`.
3. For the vertical line with zero noise and `x>0`, one step should reduce `x`; for `x<0`, one step should increase `x`.
4. For Gaussian measurement noise with a fixed seed, output is reproducible.
5. Local projection should choose endpoints if Newton fails or leaves the interval.
6. Ellipse speed formula must match `sqrt(a^2 sin^2 t + b^2 cos^2 t)`.
7. The second-order ellipse lookahead formula must use `f(t)^4` in the denominator.

If adding tests is too much for the current task, at least add a small executable or command that runs line and ellipse smoke tests and writes CSV files.

## Plotting and analysis

C++ should produce CSV data. Plotting can be done separately, preferably with Python scripts if needed.

Do not bake plotting into the core C++ simulation unless requested. If adding plotting helpers, put them under `scripts/` and keep them optional.

## Documentation expectations

When implementing formulas, add short comments explaining:

- What physical/modeling role the formula plays.
- Which parameter corresponds to the mathematical draft.
- Any approximation being used, especially for lookahead `T`.

Do not over-comment obvious C++ syntax.

## Safety and repository hygiene

- Do not commit generated binaries.
- Do not commit large output CSVs unless explicitly requested.
- Put generated experiment output under `outputs/`.
- Keep random seeds visible in output metadata.
- Do not delete or rewrite the mathematical draft unless explicitly asked.
- If modifying formulas from the draft, document the reason.

## Suggested first Codex task

If asked to start implementation from scratch, first create:

1. `Vec2` utilities.
2. `Curve2D` interface with line and ellipse implementations.
3. `SimulationConfig` and `State`.
4. `wrap_angle`, expected-direction calculation, direction/speed update.
5. Local projection with Newton candidate plus endpoint fallback.
6. CSV output.
7. A simple line experiment and an ellipse experiment.

After this foundation works, add comparative experiments for:

- first-order vs second-order ellipse lookahead,
- same-noise comparison runs,
- sparse sampling,
- measurement noise.
