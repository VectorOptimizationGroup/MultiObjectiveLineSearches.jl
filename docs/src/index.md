# MultiObjectiveLineSearches.jl

MultiObjectiveLineSearches provides line searches for multiobjective
optimization. Given a point, a search direction and the objective values at
the point, a line search looks for a step size along the direction that
satisfies an acceptance condition on all objectives.

The available method is [`BackTracking`](@ref): backtracking for an
Armijo-type sufficient decrease condition, with safeguarded quadratic
interpolation.

## Contents

- [Using a line search](@ref) explains evaluation, buffers and result handling.
- [BackTracking](@ref backtracking) presents the formulation, literature,
  safeguards and an executable example.
- [API Reference](@ref) documents the shared configuration and result types.
- [References](@ref) lists the cited publications.

## Installation

The package is not registered. From a separate Julia environment, use
`Pkg.develop(path="/path/to/MultiObjectiveLineSearches.jl")` with your checkout.
Julia 1.11 or later is required.
