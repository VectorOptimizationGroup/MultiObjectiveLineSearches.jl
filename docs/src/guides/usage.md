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
component of `F_trial`, must not modify `x_trial` or the input vectors, and
must return:

- `nothing` after a completed evaluation;
- an error message when it explicitly reports an evaluation failure.

An exception from the evaluator propagates. A `NaN` or infinite value
returned with `nothing` is not an error: the trial is rejected, and
BackTracking halves the step if its limits allow another trial. Returning
anything other than `nothing` or a string raises `ArgumentError`.

| Argument | Content |
|:---|:---|
| `x`, `d` | Current point and direction |
| `Fx` | All objectives evaluated at `x` |
| `slopes` | One slope per objective, used only to propose reduced steps |
| `factor` | Finite, negative scalar in the acceptance condition |
| `x_trial`, `F_trial` | Writable trial point and objective buffers |

All six vectors must have the configured element type. The axes of `x`, `d`
and `x_trial` must match; those of `Fx`, `slopes` and `F_trial` must also match.
`x`, `d`, `Fx` and `slopes` must be nonempty and contain only finite values. Dimension mismatches raise
`DimensionMismatch`; invalid values, element types or storage overlap raise
`ArgumentError` before any evaluation. `factor` is converted to the configured
type and must remain finite and negative.

`factor` enters the acceptance condition; `slopes` are used only to propose
reduced steps. The search computes no derivatives and does not check how these
numbers were obtained: their meaning is set by the calling method. See
[Armijo-type conditions in the literature](@ref armijo-literature) for choices
made in the cited works.

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
solution. BackTracking does not modify `x`, `d`, `Fx` or `slopes`. To retain an
accepted point in a history, copy it before the buffers are reused. [`LineSearchResult`](@ref) contains
scalar diagnostics, not copies of the objective or point vectors.

| `reason` | Meaning |
|:---|:---|
| `:success` | All objectives satisfy the acceptance condition |
| `:max_iterations` | The last trial was rejected after `max_iters` reductions |
| `:minimum_step` | A reduced step is at or below `minstep`; it was not evaluated |
| `:evaluation_failure` | The evaluator returned an error message |
| `:indistinguishable_point` | `x + alpha*d` equals `x` in floating point; not evaluated |

`alpha` is an accepted step only on success. Otherwise it is the last candidate,
possibly not evaluated. `nfev` counts actual evaluator calls, including a call
that reports an error; `iterations` counts reductions. The initial trial uses
zero reductions. The solver determines how a failed search affects its own
termination.
