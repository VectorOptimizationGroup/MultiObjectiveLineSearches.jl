# Quick Start

After [installing MultiObjectiveLineSearches.jl](@ref Installation), define an
in-place evaluator for ``F_1(x)=(x-1)^2`` and ``F_2(x)=(x-1)^2+1``:

```jldoctest quickstart
julia> using MultiObjectiveLineSearches

julia> function objectives!(values, x)
           values[1] = (x[1] - 1)^2
           values[2] = (x[1] - 1)^2 + 1
           return nothing
       end;

julia> x = [0.0]; d = [3.0]; Fx = [1.0, 2.0]; slopes = [-6.0, -6.0];
```

Prepare once using the evaluator and initial point/values:

```jldoctest quickstart
julia> search = configure_linesearch(BackTracking(), objectives!, x, Fx);

julia> result = perform_linesearch!(search, x, d, Fx;
                                  slopes, factor=maximum(slopes), maxstep=1);

julia> result.success, result.alpha
(true, 0.3333333333333333)

julia> result.x_trial, result.F_trial
([1.0], [0.0, 1.0])

julia> if result.success
           copyto!(x, result.x_trial)
           copyto!(Fx, result.F_trial)
       end;

julia> x, Fx
([1.0], [0.0, 1.0])
```

Reuse `search` with new iteration data. Copy accepted results before the next
call: their arrays borrow reusable buffers. On failure, keep the previous state.
See [Using a line search](@ref) for the full contract and
[BackTracking](@ref backtracking) for the mathematical explanation.
