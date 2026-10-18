# MultiObjectiveLineSearches.jl

Multiobjective line searches with reusable preparation and in-place execution.
The available method is [`BackTracking`](@ref), with an Armijo-type decrease
condition and safeguarded quadratic interpolation.

- [Installation](@ref) and [Quick Start](@ref): get started.
- [Using a line search](@ref): inputs, evaluators and result handling.
- [Extending the protocol](@ref): implement a method with its own state.
- [BackTracking](@ref backtracking): mathematics and literature comparisons.
- [API Reference](@ref) and [References](@ref).
