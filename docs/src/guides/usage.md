# Using a line search

A solver chooses the direction and supplies the objective values at its current
point. The line search evaluates trial points and returns a step plus a
diagnostic. See [BackTracking](@ref backtracking) for the mathematical test and
a complete executable example.

## Configure once, reuse across iterations

Create a [`BackTracking`](@ref) specification and call
[`configure_linesearch`](@ref) with the element type of the solver's vectors.
The returned configuration is reusable. [`InitialStatic`](@ref) starts every
search at its configured step; it does not remember the previous accepted step.

Allocate `x_trial = similar(x)` and `F_trial = similar(Fx)` once and reuse them.
They must not share storage with each other or with any input vector, including
through views.

## Supply the evaluation and directional data

The evaluator has the form `evaluate!(F_trial, x_trial)`. It must fill every
component of `F_trial`, leave the point and input vectors unchanged, and return:

- `nothing` after a completed evaluation;
- an error message when it explicitly reports an evaluation failure.

An exception from the evaluator propagates. A nonfinite objective returned with
`nothing` is instead a trial that BackTracking can reject and retry at a smaller
step. Returning anything other than `nothing` or a string is an API error.

| Argument | Content |
|:---|:---|
| `x`, `d` | Current point and direction |
| `Fx` | All objectives evaluated at `x` |
| `slopes` | Component slopes used by interpolation |
| `factor` | Finite, negative scalar used for sufficient decrease |
| `x_trial`, `F_trial` | Writable trial point and objective buffers |

All six vectors must have the configured element type. The axes of `x`, `d`
and `x_trial` must match; those of `Fx`, `slopes` and `F_trial` must also match.
Input vectors must be nonempty and finite. Dimension mismatches raise
`DimensionMismatch`; invalid values, element types or storage overlap raise
`ArgumentError` before any evaluation. `factor` is converted to the configured
type and must remain finite and negative.

`slopes` and `factor` have different roles. For smooth objectives a common
choice is `slopes = JF*d`, `factor = maximum(slopes)`. For composite objectives,
use the direction and gap appropriate to the solver, as explained in
[BackTracking](@ref backtracking).

## Consume the result

If `result.success` is true, the two trial buffers hold the accepted point and
objective values. The caller can copy them into its current state:

```julia
if result.success
    copyto!(x, x_trial)
    copyto!(Fx, F_trial)
else
    # Keep x and Fx, and handle result.reason in the calling solver.
end
```

On failure, trial buffers are scratch data and do not represent a new accepted
solution. Inputs remain unchanged. To retain an accepted point in a history,
copy it before the buffers are reused. [`LineSearchResult`](@ref) contains
scalar diagnostics, not copies of the objective or point vectors.

| `reason` | Meaning |
|:---|:---|
| `:success` | All components passed the acceptance test |
| `:max_iterations` | No accepted step within the reduction budget |
| `:minimum_step` | A reduced candidate reached or crossed `minstep` |
| `:evaluation_failure` | The evaluator returned an error message |
| `:indistinguishable_point` | The trial point rounded to the current point |

`alpha` is an accepted step only on success. Otherwise it is the last candidate,
possibly not evaluated. `nfev` counts actual evaluator calls, including a call
that reports an error; `iterations` counts reductions. The initial trial uses
zero reductions. The solver determines how a failed search affects its own
termination.
