# PhysicalParticles.jl

Physical particle types for scientific simulation

Referred to [GeometicalPredicates](https://github.com/JuliaGeometry/GeometricalPredicates.jl)

## Installation

```julia
]add PhysicalParticles
```

or

```julia
using Pkg; Pkg.add("PhysicalParticles")
```

## Usage

### Mathematical Points



### Physical vectors

```julia
julia> p1 = PosAstro(1.0,2.0,3.0)
PhysicalVector3D(1.0 kpc, 2.0 kpc, 3.0 kpc)

julia> v = VelAstro(5.0,6.0,7.0)
PhysicalVector3D(5.0 kpc Gyr^-1, 6.0 kpc Gyr^-1, 7.0 kpc Gyr^-1)

julia> p2 = p1 + v * 5u"s"
PhysicalVector3D(3.0856775814913696e19 m, 6.1713551629827375e19 m, 9.257032744474105e19 m)
```
