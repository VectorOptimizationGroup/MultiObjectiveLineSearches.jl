# Extending the protocol

New methods implement the same two functions with method-specific prepared
state; no central dispatcher or common prepared supertype is required.

1. Define a specification subtype of [`AbstractLineSearchMethod`](@ref).
2. Implement `configure_linesearch(spec::YourMethod, evaluate!, x, Fx; ...)`:
   convert/validate fixed options, check required callback applicability without
   evaluation, and allocate independent buffers and any history seed.
3. Implement `perform_linesearch!(search::YourPrepared, x, d, Fx;
   maxstep=Inf, method_data...)`: validate new inputs/aliasing and bounds,
   apply the acceptance rule, and return [`LineSearchResult`](@ref).
   Reject unsupported keywords.

Follow [Using a line search](@ref) for buffer ownership and history commits.
Store only what the method needs; BackTracking needs neither derivatives nor history.

## Derivative callbacks

Supply optional capabilities during preparation:

```julia
configure_linesearch(spec, evaluate!, x, Fx;
                     derivative! = nothing, value_derivative! = nothing)
```

`derivative!(D, y, d)` fills per-objective directional derivatives with the
objective axes/type; `value_derivative!(F, D, y, d)` fills both values and
derivatives. Their return/error contract is the same as `evaluate!`.
Keep a separate value evaluator even when a joint callback exists.

Methods requiring derivatives reject preparation if neither capability exists.
Check applicability again if new direction containers can affect it.
BackTracking ignores the optional callbacks. Each new method must specify its
acceptance rule; derivative support alone does not establish suitability for
nonsmooth/composite objectives.

| Request | Callback | Counter increments |
| --- | --- | --- |
| Values | `evaluate!` | `nfev` |
| Derivatives, separate callback available | `derivative!` | `ndev` |
| Derivatives, only joint available | joint, with value scratch | `nfev`, `ndev` |
| Both, joint available | joint | `nfev`, `ndev` |
| Both, no joint | values, then derivatives | each requested counter; stop if values fail or are nonfinite |

Increment counters before calling, including reported failures. A joint call
adds one to both `nfev` and `ndev`, even if one output is discarded. Reused or
supplied data do not count again; internal callback work/cache is not tracked.
The sum is not a count of callback calls or computational cost. Nonfinite output
follows the method's acceptance/failure rule.

## Example extensions

The checkout's `test/protocol_fixtures.jl` contains two test-only examples:
`DerivativeProbe` accepts finite values no greater than the base and nonpositive
finite directional derivatives; `HistoryProbe` accepts finite values no greater
than their componentwise history maximum. They illustrate the protocol, not
production Wolfe or nonmonotone algorithms.

```jldoctest extensions
julia> using MultiObjectiveLineSearches

julia> include(joinpath(pkgdir(MultiObjectiveLineSearches), "test", "protocol_fixtures.jl"));

julia> using .ProtocolFixtures: DerivativeProbe, HistoryProbe

julia> values!(F, y) = (fill!(F, (y[1] - 1)^2); nothing);

julia> directional!(D, y, d) = (fill!(D, 2 * (y[1] - 1) * d[1]); nothing);

julia> function joint!(F, D, y, d)
           values!(F, y)
           return directional!(D, y, d)
       end;

julia> function accept_steps!(search, x, Fx, steps; maxstep=1)
           last_result = nothing
           for (d, data) in steps
               result = perform_linesearch!(search, x, d, Fx; maxstep, data...)
               result.success || return result
               copyto!(x, result.x_trial)
               copyto!(Fx, result.F_trial)
               last_result = result
           end
           return last_result
       end;
```

The same consumer handles all methods; callers supply method-specific data:

```jldoctest extensions
julia> plain_steps = [([0.5], (;)), ([0.25], (;))];

julia> bt_steps = [([0.5], (slopes=[-1.0, -1.0], factor=-1.0)),
                  ([0.25], (slopes=[-0.25, -0.25], factor=-0.25))];

julia> function demo(spec, steps; capabilities...)
           x = [0.0]; Fx = [1.0, 1.0]
           search = configure_linesearch(spec, values!, x, Fx; capabilities...)
           result = accept_steps!(search, x, Fx, steps)
           return x, Fx, result, search
       end;

julia> xb, Fb, rb, _ = demo(BackTracking(), bt_steps);

julia> xs, Fs, rs, _ = demo(DerivativeProbe(), plain_steps; derivative! = directional!);

julia> xj, Fj, rj, _ = demo(DerivativeProbe(), plain_steps; value_derivative! = joint!);

julia> xh, Fh, rh, history = demo(HistoryProbe(), plain_steps);

julia> all(x -> x == [0.75], (xb, xs, xj, xh)), all(F -> F == [0.0625, 0.0625], (Fb, Fs, Fj, Fh))
(true, true)

julia> (rb.nfev, rb.ndev), (rs.nfev, rs.ndev), (rj.nfev, rj.ndev)
((1, 0), (1, 2), (2, 2))

julia> history.accepted, history.used
(2, 3)

julia> try
           configure_linesearch(DerivativeProbe(), values!, [0.0], [1.0, 1.0])
       catch error
           error isa ArgumentError
       end
true
```

Counts refer to the last invocation. The history seed occupies one slot;
two accepted searches add two commits. Behavioral tests are in `test/protocol.jl`.
