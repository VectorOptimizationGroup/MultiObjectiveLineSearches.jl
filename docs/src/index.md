# MultiObjectiveLineSearches.jl

MultiObjectiveLineSearches provides line searches for vector-valued objectives.
Given a point, a direction and objective values, a search chooses a step that
satisfies sufficient decrease in every objective.

The available method is [`BackTracking`](@ref), an Armijo search with
safeguarded quadratic interpolation. It supports smooth and composite
objectives through a user-supplied decrease factor.

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
