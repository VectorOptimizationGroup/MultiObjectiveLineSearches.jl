# [BackTracking](@id backtracking)

`BackTracking` searches along a given direction for a step size that satisfies
an Armijo-type sufficient decrease condition in every objective. Rejected steps
are reduced by safeguarded quadratic interpolation.

## Overview

Let ``F=(F_1,\ldots,F_m):\mathbb{R}^n\to\mathbb{R}^m`` be the objective map,
let ``x`` be the current point and let ``d(x)`` be a search direction at ``x``.
Along ``d(x)``, write

```math
\phi_i(\alpha)=F_i(x+\alpha d(x)),\qquad i=1,\ldots,m.
```

`BackTracking` uses only the values ``\phi_i(\alpha)`` returned by the caller's
evaluator at trial steps, together with two inputs from the caller: a scalar
decrease factor ``\eta<0`` (`factor`) and one slope ``s_i`` per objective
(`slopes`). It computes no derivatives and makes no assumption about ``F``,
such as differentiability or convexity. The meaning of ``\eta`` and ``s_i`` is
set by the method that calls the search;
[Armijo-type conditions in the literature](@ref armijo-literature) shows the
choices made in the cited works.

| Option | Default | Meaning |
|:---|:---|:---|
| `initial_step` | `InitialStatic(alpha=1.0)` | First trial step of each search |
| `c_1` | `1e-4` | Sufficient decrease coefficient |
| `tau_1` | `0.05` | Lower safeguard ``\tau_1`` of the quadratic interpolation |
| `tau_2` | `0.95` | Upper safeguard ``\tau_2`` of the quadratic interpolation |
| `minstep` | `1e-15` | A reduced step at or below this value ends the search |
| `maxstep` | `1.0` | Upper bound for the initial step |
| `max_iters` | `1000` | Maximum number of step reductions |
| `verbose` | `0` | `1` logs the final diagnostic with `@info` |

[`configure_linesearch`](@ref) converts the options to the element type `T` and
checks that ``0<c_1<1``, that ``0<\tau_1\leq1/2\leq\tau_2<1``, that
`minstep`, `maxstep` and the initial step are finite, that
`0 < minstep < maxstep`, and that `minstep <= initial_step.alpha <= maxstep`.
Steps only decrease during a search, so no trial step exceeds `maxstep`.

## [Armijo-type condition](@id armijo-condition)

`BackTracking` accepts ``\alpha>0`` when every ``F_i(x+\alpha d(x))`` is finite
and

```math
F_i(x+\alpha d(x))\leq F_i(x)+c_1\alpha\eta,\qquad i=1,\ldots,m.
```

The same term ``c_1\alpha\eta`` is used for all objectives. For ``m=1`` and
``\eta=\nabla F_1(x)^\mathsf{T}d(x)``, this is the Armijo condition of scalar
optimization; see Nocedal and Wright [NW2006](@cite), equation (3.4).

The multiobjective works cited here name their extensions of this condition
differently: Fliege and Svaiter [FS2000](@cite) call it an Armijo-like rule,
Lucambio Pérez and Prudente [LPP2018](@cite) the Armijo condition for vector
optimization, and Assunção, Ferreira and Prudente [AFP2025](@cite) an
Armijo-type step size. This documentation uses *Armijo-type condition*.

When ``F`` is differentiable, the generalized directional derivative of ``F``
at ``x`` along ``d(x)`` is

```math
\mathcal{D}(x,d(x))=\max_{1\leq j\leq m}\nabla F_j(x)^\mathsf{T}d(x),
```

see [LPP2018](@cite). If ``\mathcal{D}(x,d(x))<0``, then ``d(x)`` is a
descent direction for every objective.

`BackTracking` does not verify that the condition can be satisfied: with an
arbitrary negative `factor`, the search may end with `:minimum_step` or
`:max_iterations`. For each condition in
[Armijo-type conditions in the literature](@ref armijo-literature), the cited
work shows that, under its assumptions, all sufficiently small steps are
accepted.

## [Armijo-type conditions in the literature](@id armijo-literature)

The works below differ in the decrease term on the right-hand side of the
condition. In all of them, ``c_1\in(0,1)``.

### One slope per objective

Fliege and Svaiter [FS2000](@cite), Section 4, assume that ``F`` is
differentiable and that ``\nabla F_i(x)^\mathsf{T}d(x)<0`` for all ``i``, and
accept ``\alpha`` when

```math
F_i(x+\alpha d(x))\leq F_i(x)+c_1\alpha\nabla F_i(x)^\mathsf{T}d(x),
\qquad i=1,\ldots,m,
```

their equation (7). The search starts at ``\alpha=1`` and halves ``\alpha``
until (7) holds. Their Lemma 4 shows that (7) holds strictly for all
sufficiently small ``\alpha>0``.

Each objective has its own slope ``s_i=\nabla F_i(x)^\mathsf{T}d(x)`` on the
right-hand side. `BackTracking` does not implement this rule, since it uses one
scalar ``\eta`` for all objectives. With ``\eta=\mathcal{D}(x,d(x))`` we have
``s_i\leq\eta``, so (7) implies the `BackTracking` condition; the converse does
not hold in general.

### The generalized directional derivative

Lucambio Pérez and Prudente [LPP2018](@cite) work in vector optimization, where
vectors are compared through the partial order induced by a cone. In the
multiobjective case this order is componentwise, and their Armijo condition,
equation (6), reduces to equation (7) of Fliege and Svaiter [FS2000](@cite).
Their Definition 3.1 introduces Wolfe conditions for vector optimization. With
the choices the authors indicate for multiobjective problems, the sufficient
decrease part of these conditions, (7a) and (8a), uses the generalized
directional derivative as the common term:

```math
F_i(x+\alpha d(x))\leq F_i(x)+c_1\alpha\mathcal{D}(x,d(x)),
\qquad i=1,\ldots,m.
```

This is the `BackTracking` condition with ``\eta=\mathcal{D}(x,d(x))``.
The proof of their Proposition 3.2 shows that it holds for all sufficiently
small ``\alpha>0`` when ``F`` is continuously differentiable,
``\mathcal{D}(x,d(x))<0`` and ``F`` is bounded below along ``d(x)``.
`BackTracking` checks only this decrease condition, not the curvature
conditions (7b) or (8b).

### Composite objectives

Assunção, Ferreira and Prudente [AFP2025](@cite) consider objectives
``F_i=f_i+h_i``, where ``f_i`` is continuously differentiable and ``h_i`` is
proper, convex and lower semicontinuous. At an iterate ``x``, their generalized
conditional gradient method computes, in equations (14) and (15),

```math
\begin{aligned}
p(x)&\in\operatorname*{argmin}_{u}\max_{1\leq j\leq m}
\{h_j(u)-h_j(x)+\nabla f_j(x)^\mathsf{T}(u-x)\},\\
\theta(x)&=\max_{1\leq j\leq m}
\{h_j(p(x))-h_j(x)+\nabla f_j(x)^\mathsf{T}(p(x)-x)\},
\end{aligned}
```

and uses the direction ``d(x)=p(x)-x``. By their Lemma 3.1,
``\theta(x)\leq0``, with equality if and only if ``x`` is Pareto critical.
Their Armijo-type step size (Section 4, steps LS0–LS2) starts at ``\alpha=1``,
accepts ``\alpha`` when

```math
F_i(x+\alpha d(x))\leq F_i(x)-c_1\alpha|\theta(x)|,
\qquad i=1,\ldots,m,
```

and otherwise chooses the next step in ``[\tau_1\alpha,\tau_2\alpha]``, with
``0<\tau_1<\tau_2<1``. This is the `BackTracking` condition with
``\eta=-|\theta(x)|``. Their Proposition 4.1 shows that, when ``\theta(x)<0``,
the condition holds for all ``\alpha`` in an interval ``(0,\bar\alpha]`` with
``\bar\alpha\leq1``. With the default options, every trial step lies in
``(0,1]``, the range of step sizes in their Algorithm 1. Every reduction made
by `BackTracking` lies in ``[\tau_1\alpha,\tau_2\alpha]``; see
[Quadratic interpolation and safeguards](@ref quadratic-interpolation).
The defaults ``c_1=10^{-4}``, ``\tau_1=0.05`` and ``\tau_2=0.95`` are the
values reported in their Section 5, where the line search is coded by
quadratic interpolation of the coordinate functions.

The objectives ``F_i`` need not be differentiable here, so the slopes used for
interpolation are chosen by the caller, for example the directional
derivatives ``\nabla f_i(x)^\mathsf{T}d(x)`` of the differentiable parts.

## [Quadratic interpolation and safeguards](@id quadratic-interpolation)

Suppose a trial with finite values is rejected, and let ``j`` be the first
index, in increasing order, whose condition fails. With the current step
``\alpha`` and the caller's slope ``s_j``, `BackTracking` forms

```math
q_j(t)=\phi_j(0)+s_jt+b_jt^2,
\qquad
b_j=\frac{\phi_j(\alpha)-\phi_j(0)-\alpha s_j}{\alpha^2},
```

so that ``q_j(0)=\phi_j(0)``, ``q_j'(0)=s_j`` and ``q_j(\alpha)=\phi_j(\alpha)``.
The stationary point of ``q_j`` is

```math
\widehat\alpha
=-\frac{s_j\alpha^2}
        {2\,[\phi_j(\alpha)-\phi_j(0)-\alpha s_j]}.
```

The new step is ``\widehat\alpha`` if ``s_j<0``, the denominator is finite and
nonzero, ``\widehat\alpha`` is finite, and

```math
\tau_1\alpha\leq\widehat\alpha\leq\tau_2\alpha,
```

where ``\tau_1`` and ``\tau_2`` are the options `tau_1` and `tau_2`.

In that case ``\widehat\alpha>0`` and ``s_j<0`` imply ``b_j>0``, so
``\widehat\alpha`` is the minimizer of ``q_j``. Otherwise the new step is
``\alpha/2``: a proposal outside the interval is replaced, not clamped to it.
If some trial value is `NaN` or infinite, the step is halved without
interpolation. Since ``\tau_1\leq1/2\leq\tau_2``, halving also lies in the
interval, so every reduction gives a new step in ``[\tau_1\alpha,\tau_2\alpha]``;
with the defaults, between 5% and 95% of the previous one. With
`tau_1 = tau_2 = 0.5`, every reduction halves the step, as in the search of
Fliege and Svaiter [FS2000](@cite).

When ``s_j=\phi_j'(0)``, ``q_j`` and ``\widehat\alpha`` are the quadratic
interpolant and its minimizer in Nocedal and Wright [NW2006](@cite),
Section 3.5, equations (3.57) and (3.58). Otherwise ``q_j`` still matches the two
values and the slope ``s_j``. The slopes only affect the proposed steps; the
acceptance test uses the objective values and ``\eta``.

Nocedal and Wright continue with cubic interpolation after the first quadratic
trial, and reset a trial to half its predecessor when it is too close to or
too much smaller than that predecessor. `BackTracking` uses the quadratic
model at every reduction of a trial with finite values, with the interval
above. Choosing the first violating objective and the fallback rules above are
choices of this implementation; the cited works do not prescribe them.

## Trial sequence and stopping

Starting from `alpha = initial_step.alpha`:

1. Form ``x+\alpha d(x)``. If it equals ``x`` in floating-point arithmetic,
   return `:indistinguishable_point` without evaluating it.
2. Evaluate the objectives. If the evaluator reports an error, return
   `:evaluation_failure`.
3. If all values are finite and satisfy the Armijo-type condition, return
   `:success`.
4. If `max_iters` reductions have been made, return `:max_iterations`.
5. Reduce the step: halve it if some value is not finite, otherwise apply the
   safeguarded quadratic proposal.
6. If the new step is at or below `minstep`, return `:minimum_step` without
   evaluating it. Otherwise go to step 1.

With `max_iters=0`, only the initial step is evaluated. An initial step equal
to `minstep` is evaluated, since the threshold applies to reduced steps.
The limits `max_iters` and `minstep`, the rounding check and the handling of
nonfinite values belong to this implementation; they are not part of the line
searches of Fliege and Svaiter [FS2000](@cite) or of Assunção, Ferreira and
Prudente [AFP2025](@cite).

## Usage

Consider ``F=(F_1,F_2):\mathbb{R}\to\mathbb{R}^2`` with

```math
F_1(x)=(x-1)^2,\qquad F_2(x)=(x-1)^2+1,
```

the current point ``x=0`` and the direction ``d(x)=3``. Then ``F(0)=(1,2)``
and ``\nabla F_1(0)=\nabla F_2(0)=-2``, so the slopes are
``s_i=\nabla F_i(x)^\mathsf{T}d(x)=-6``. The caller passes
``\eta=\mathcal{D}(x,d(x))=-6``, the choice of Lucambio Pérez and Prudente
[LPP2018](@cite).

The first trial, ``\alpha=1``, gives ``F(3)=(4,5)`` and is rejected by
``F_1``. Since ``F_1`` is quadratic and ``s_1=\phi_1'(0)``, the model ``q_1``
equals ``\phi_1``, and its minimizer ``\widehat\alpha=1/3`` is accepted with
``F(1)=(0,1)``. In the code, `objectives!` fills `values` with ``F(x)``.

```jldoctest backtracking_usage
julia> using MultiObjectiveLineSearches

julia> function objectives!(values, x)
           values[1] = (x[1] - 1)^2
           values[2] = (x[1] - 1)^2 + 1
           return nothing
       end;

julia> x = [0.0]; d = [3.0]; Fx = [1.0, 2.0]; slopes = [-6.0, -6.0];

julia> search = configure_linesearch(BackTracking(), Float64);

julia> x_trial = similar(x); F_trial = similar(Fx);

julia> result = perform_linesearch!(search, objectives!, x, d, Fx, slopes,
                                  x_trial, F_trial; factor=maximum(slopes));

julia> (result.success, result.nfev, result.iterations, result.alpha ≈ 1/3)
(true, 2, 1, true)

julia> (x_trial ≈ [1.0], F_trial ≈ [0.0, 1.0], x == [0.0])
(true, true, true)
```

For the evaluator contract and buffer ownership, see [Using a line search](@ref).

## Method reference

```@docs
MultiObjectiveLineSearches.BackTracking
MultiObjectiveLineSearches.perform_linesearch!
```
