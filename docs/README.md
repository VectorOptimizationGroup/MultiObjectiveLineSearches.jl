# Building the documentation

From the repository root:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The build checks exported docstrings and runs doctests, then writes
`docs/build/index.html`. Extension examples use `test/protocol_fixtures.jl`
from this checkout. Local builds skip deployment.

Use the same Julia executable for setup and build. When testing another minor
version, resolve a separate docs environment instead of reusing its manifest.

To serve the generated site with LiveServer installed:

```bash
julia -e 'using LiveServer; serve(dir="docs/build")'
```

Open `http://localhost:8000/`.
