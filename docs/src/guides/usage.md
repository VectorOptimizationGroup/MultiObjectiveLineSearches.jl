# Using a line search

Prepare once, then reuse the search with new iteration data:

```julia
search = configure_linesearch(spec, evaluate!, x, Fx)
result = perform_linesearch!(search, x, d, Fx; maxstep=Inf, method_data...)
```

See [Quick Start](@ref) for a complete example and
[Extending the protocol](@ref) for derivative callbacks and new methods.

## Preparation and inputs

Preparation converts fixed options to `T = eltype(x)` and allocates independent
work buffers without calling the evaluator. `T` must be a concrete floating
type, `eltype(Fx) === T`, and both prototypes must be finite and nonempty.
Mixed types raise `ArgumentError`; there is no automatic promotion.

Reuse with matching axes/types. Changes to dimensions, type, problem or callback
capabilities require a new preparation. Do not modify prepared storage or share
one prepared object across concurrent calls. Captured callback state remains
the caller's responsibility; keep BigFloat precision/rounding fixed during a run.

BackTracking requires these iteration inputs:

| Input | Requirement |
| --- | --- |
| `x`, `d` | Finite point and direction, with prepared point axes/type |
| `Fx` | Finite objectives at `x`, with prepared value axes/type |
| `slopes` | Finite vector with objective axes/type; used only for interpolation |
| `factor` | Real scalar, finite and negative after conversion; controls acceptance |
| `maxstep` | Nonnegative real bound; `Inf` adds no restriction |

New data are validated each call; fixed options and buffers are reused.
Inputs may alias each other but must not alias prepared outputs, including
through views. Invalid axes raise `DimensionMismatch`; invalid types, values
or aliasing raise `ArgumentError`. The search trusts `Fx` and the mathematical
meaning of `slopes`/`factor`; it computes no derivatives.

## Evaluator

`evaluate!(F, y)` overwrites all objective values without modifying `y` or
captured inputs. Return `nothing` on completion or an `AbstractString` on
failure; other returns raise `ArgumentError`, and exceptions propagate.

Nonfinite values returned with `nothing` cause BackTracking to reject and halve
the step if limits allow. Trial-point overflow can reach the evaluator: restrict
the domain or report failure if the evaluator requires finite points.

## Step limit

The per-call `maxstep` is converted without rounding above the requested bound.
Zero is valid; negative, NaN or non-real bounds raise `ArgumentError`.
BackTracking starts at the smaller of its initial step and this bound.
A bound below `minstep` returns `:step_limit` without evaluation; equality
permits the initial trial. See [Trial sequence and stopping](@ref) for subsequent
reductions and termination reasons.

## Results and failure

On success, advance the consumer using `result.x_trial` and `result.F_trial`.
These arrays borrow reusable buffers: copy them before the next call, and never
use them as inputs or modify them through a retained result. Scalar diagnostics
remain valid; retained arrays are not snapshots.

```jldoctest usage
julia> using MultiObjectiveLineSearches

julia> function objectives!(values, x)
           values[1] = (x[1] - 1)^2
           values[2] = (x[1] - 1)^2 + 1
           return nothing
       end;

julia> x = [0.0]; d = [3.0]; Fx = [1.0, 2.0]; slopes = [-6.0, -6.0];

julia> search = configure_linesearch(BackTracking(), objectives!, x, Fx);

julia> result = perform_linesearch!(search, x, d, Fx; slopes, factor=-6.0);

julia> if result.success
           copyto!(x, result.x_trial)
           copyto!(Fx, result.F_trial)
       end;

julia> x, Fx
([1.0], [0.0, 1.0])

julia> bounded = perform_linesearch!(search, x, d, Fx;
                                   slopes, factor=-6.0, maxstep=0);

julia> bounded.success, bounded.reason, bounded.alpha, bounded.nfev
(false, :step_limit, 0.0, 0)
```

Failure leaves the consumer's accepted state unchanged. Result arrays are
scratch and `alpha` is a diagnostic candidate, possibly unevaluated or zero.
The consumer decides how to handle termination; [`LineSearchResult`](@ref)
describes its fields and counters. BackTracking has `ndev == 0`.

To change the evaluator, prepare again:

```jldoctest usage
julia> failing!(values, y) = "objective undefined at the trial point";

julia> failed_search = configure_linesearch(BackTracking(), failing!, x, Fx);

julia> failed = perform_linesearch!(failed_search, x, d, Fx; slopes, factor=-6.0);

julia> failed.success, failed.reason, failed.nfev
(false, :evaluation_failure, 1)

julia> failed.message
"objective undefined at the trial point"

julia> x, Fx
([1.0], [0.0, 1.0])
```

## Acceptance and history

Success commits acceptance. If the consumer vetoes it, discard that prepared
execution before continuing. Methods with history copy the initial seed, require
each call to start at the last committed point/values, and commit once per
success. Rejected trials, failures and exceptions do not commit history;
callback side effects and scratch storage need not roll back.
BackTracking carries no history and performs no continuity checks.
