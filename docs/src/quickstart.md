# Quick Start

After [installing MultiObjectiveLineSearches.jl](@ref Installation), load the
package and configure a line search for vectors of `Float64`:

```jldoctest quickstart
julia> using MultiObjectiveLineSearches

julia> search = configure_linesearch(BackTracking(), Float64);
```

The configured search is reusable across the iterations of a solver.

The objectives are evaluated by an in-place function that fills its first
argument and returns `nothing`. This example uses
``F_1(x)=(x-1)^2`` and ``F_2(x)=(x-1)^2+1``:

```jldoctest quickstart
julia> function objectives!(values, x)
           values[1] = (x[1] - 1)^2
           values[2] = (x[1] - 1)^2 + 1
           return nothing
       end;
```

Supply the current point, the search direction, the objective values at the
point, one slope per objective and the decrease factor, together with buffers
for the trial point and its objective values:

```jldoctest quickstart
julia> x = [0.0]; d = [3.0]; Fx = [1.0, 2.0]; slopes = [-6.0, -6.0];

julia> x_trial = similar(x); F_trial = similar(Fx);

julia> result = perform_linesearch!(search, objectives!, x, d, Fx, slopes,
                                  x_trial, F_trial; factor = maximum(slopes));

julia> result.success, result.alpha
(true, 0.3333333333333333)

julia> x_trial, F_trial
([1.0], [0.0, 1.0])
```

On success, `x_trial` and `F_trial` hold the accepted point and its objective
values; `x`, `d`, `Fx` and `slopes` are not modified.

Continue with [Using a line search](@ref) for the evaluator, buffer and result
contracts, or [BackTracking](@ref backtracking) for the Armijo-type condition,
its sources and the choice of `slopes` and `factor`.
