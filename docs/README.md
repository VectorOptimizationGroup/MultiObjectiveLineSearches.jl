# MultiObjectiveLineSearches.jl Documentation

The Documenter.jl source lives in `docs/src/`.

To build the local HTML documentation from the repository root:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The generated site is written to `docs/build/index.html`.

To browse it through a local server:

```bash
julia -e 'using LiveServer; serve(dir="docs/build")'
```

Then open `http://localhost:8000/`.

The build runs every `jldoctest` block and fails if an example errors or its
output differs from the documented one.
