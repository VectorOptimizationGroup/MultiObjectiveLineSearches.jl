# [BackTracking](@id backtracking)

`BackTracking` finds a positive step satisfying componentwise Armijo decrease.
This page relates that condition to the smooth formulations of Fliege and
Svaiter [FS2000](@cite) and Lucambio Pérez and Prudente [LPP2018](@cite), and
the composite formulation of Assunção, Ferreira and Prudente [AFP2025](@cite).
The interpolation model follows the scalar construction in Nocedal and Wright,
Section 3.5 [NW2006](@cite).

## Overview

Let ``F:\mathbb{R}^n\to\mathbb{R}^m`` denote the objective vector, let ``x``
be the current point and let ``d`` be a search direction. Define

```math
\phi_i(\alpha)=F_i(x+\alpha d),\qquad i=1,\ldots,m.
```

The search changes only ``\alpha``. The caller supplies ``d``, the objective
values ``F(x)``, interpolation slopes, and a scalar decrease factor.

| Option | Default | Meaning |
|:---|:---|:---|
| `initial_step` | `InitialStatic(alpha=1.0)` | First trial of each search |
| `c_1` | `1e-4` | Sufficient-decrease coefficient, ``0<c_1<1`` |
| `minstep` | `1e-15` | Lower stopping threshold for a reduced step |
| `maxstep` | `1.0` | Upper bound on the initial step |
| `max_iters` | `1000` | Maximum number of reductions |
| `verbose` | `0` | `1` prints the final diagnostic |

All steps must be finite and positive, with
`minstep < maxstep` and `minstep <= initial_step.alpha <= maxstep`.
The quadratic safeguards are fixed at ``\tau_1=0.05`` and ``\tau_2=0.95``.

## Sufficient decrease and its sources

### Componentwise directional derivatives

For differentiable ``F``, write ``s_i=\nabla F_i(x)^\mathsf{T}d``.
Fliege and Svaiter use the componentwise condition

```math
F_i(x+\alpha d)\leq F_i(x)+c_1\alpha s_i,
\qquad i=1,\ldots,m.
```

Their Section 4, equation (7), starts from step one and repeatedly halves it.
When all ``s_i<0``, differentiability ensures that sufficiently small positive
steps satisfy the condition (Lemma 4) [FS2000](@cite).

Here the right-hand side has a **different slope for each objective**.
`BackTracking` instead takes one scalar `factor`; passing `slopes` does not
select this componentwise-slope acceptance rule.

### A common scalar factor

The implemented test is

```math
F_i(x+\alpha d)\leq F_i(x)+c_1\alpha\eta,
\qquad i=1,\ldots,m,\qquad \eta=\mathtt{factor}<0.
```

For smooth objectives, one choice is

```math
\eta=\max_{1\leq i\leq m}\nabla F_i(x)^\mathsf{T}d.
```

This is the sufficient-decrease part of the Wolfe conditions in Lucambio Pérez
and Prudente, Definition 3.1, equations (7a) and (8a), specialized to the
nonnegative orthant with the canonical basis and ``e=(1,\ldots,1)^\mathsf{T}``
[LPP2018](@cite). `BackTracking` checks decrease only; it does not impose the
curvature conditions (7b) or (8b).

For a common descent direction, every ``s_i\leq\eta<0``. Consequently the
componentwise-slope rule above implies the common-factor test, but the converse
need not hold. The common-factor test also holds for sufficiently small steps:
``s_i\leq\eta<c_1\eta`` gives strict first-order decrease in each component.
This argument depends on that relationship between ``\eta`` and the direction;
an arbitrary negative `factor` alone does not guarantee acceptance.

### Composite objectives and a subproblem gap

Write ``F_i=f_i+h_i``, with ``f_i`` differentiable and ``h_i`` convex, possibly
nonsmooth. In the generalized conditional-gradient construction, the caller
computes a point ``p`` and a gap

```math
p\in\operatorname*{argmin}_u\max_i
\{h_i(u)-h_i(x)+\nabla f_i(x)^\mathsf{T}(u-x)\},
\qquad
\theta(x)=\max_i\{h_i(p)-h_i(x)+\nabla f_i(x)^\mathsf{T}(p-x)\}.
```

For a nonstationary iterate, ``\theta(x)<0``. With ``d=p-x``, the Armijo rule is

```math
F_i(x+\alpha d)\leq F_i(x)-c_1\alpha|\theta(x)|,
\qquad 0<\alpha\leq1.
```

Thus the caller passes `factor=-abs(theta)`. Assunção, Ferreira and Prudente,
Section 4, steps LS0–LS2, start at one and require each reduction to lie in
``[\omega_1\alpha,\omega_2\alpha]`` [AFP2025](@cite). The bound ``\alpha\leq1``
keeps the trial on the segment between ``x`` and ``p``. The package's default
`maxstep=1.0` is suitable for this construction; a caller that changes the
bound must still enforce its problem's feasibility requirements.

The paper denotes the nonsmooth and smooth terms by ``g_i`` and ``h_i``;
the notation above uses ``h_i`` and ``f_i`` respectively. In a composite
application, `slopes` can contain the smooth-part values
``\nabla f_i(x)^\mathsf{T}d`` used for interpolation. They are distinct from
the subproblem gap and need not be derivatives of the complete ``F_i``.
The search neither solves the subproblem nor checks its optimality conditions.

## [Quadratic interpolation and safeguards](@id quadratic-interpolation)

Suppose a finite trial fails sufficient decrease. Choose its **first violating
component** ``j``. For the current step ``\alpha`` and supplied slope ``s_j``,
fit

```math
q_j(t)=\phi_j(0)+s_jt+b_jt^2,
\qquad
b_j=\frac{\phi_j(\alpha)-\phi_j(0)-\alpha s_j}{\alpha^2}.
```

The model matches the two values and the supplied slope:
``q_j(0)=\phi_j(0)``, ``q_j'(0)=s_j`` and
``q_j(\alpha)=\phi_j(\alpha)``. Its stationary point is

```math
\widehat\alpha
=-\frac{s_j\alpha^2}
        {2\,[\phi_j(\alpha)-\phi_j(0)-\alpha s_j]}.
```

When ``s_j`` is the actual directional derivative, these are the quadratic
interpolant and candidate of Nocedal and Wright, equations (3.57)–(3.58)
[NW2006](@cite). In the composite case, using the smooth-part slope gives an
interpolation model rather than an exact derivative match for the full objective.
The model proposes a step; the full objective vector always decides acceptance.

The implementation uses ``\widehat\alpha`` only if ``s_j<0``, the denominator
is finite and nonzero, the candidate is finite, and

```math
\tau_1\alpha\leq\widehat\alpha\leq\tau_2\alpha,
\qquad \tau_1=0.05,\quad\tau_2=0.95.
```

Otherwise it uses ``\alpha/2``. This rejects an unsuitable proposal; it does
not clamp it to an endpoint. The bounds avoid both negligible reductions and
excessively short steps relative to the current trial. Halving lies inside
this interval, so every reduction has the same safeguards.

Section 5 of Assunção, Ferreira and Prudente reports ``c_1=10^{-4}``,
``\omega_1=0.05`` and ``\omega_2=0.95``, with quadratic interpolation of coordinate
functions [AFP2025](@cite). Nocedal and Wright also describe halving an
unsuitable interpolated step, but their procedure continues with cubic models
after the first quadratic trial [NW2006](@cite). Here every finite reduction
uses a quadratic model. Selecting the first violating component and the exact
fallback checks above describe this implementation; the cited sources do not
specify all these details as one algorithm.

## Trial sequence and stopping

1. Evaluate the initial trial ``x+\alpha d``. If it rounds to ``x``, stop
   without evaluating the objective.
2. Accept only if every objective is finite and satisfies the common-factor test.
3. If the reduction budget is exhausted, return `:max_iterations`.
4. If any trial objective is `NaN` or infinite, halve the step; otherwise use
   the safeguarded quadratic proposal for the first violating component.
5. Count the reduction. If the new step is at or below `minstep`, return
   `:minimum_step` before evaluating it; otherwise repeat.

`max_iters=0` still permits the initial evaluation. An initial step equal to
`minstep` is evaluated; the threshold applies to subsequent reductions.
Explicit evaluator errors terminate immediately. Nonfinite trial values,
rounding checks and finite budgets are numerical safeguards of the package,
not additional convergence results attributed to the papers.

## Usage

This smooth example needs one quadratic reduction. The slopes are
``JF(x)d=(-6,-6)`` and the common factor is their maximum.

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

For evaluator failures and buffer ownership, see [Using a line search](@ref).

## Method reference

```@docs
MultiObjectiveLineSearches.BackTracking
MultiObjectiveLineSearches.perform_linesearch!
```
