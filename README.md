# PhysicalParticles.jl

Physical particle types for scientific simulation

Referred to [GeometicalPredicates.jl](https://github.com/JuliaGeometry/GeometricalPredicates.jl)

## Installation

```julia
]add PhysicalParticles
```

or

```julia
using Pkg; Pkg.add("PhysicalParticles")
```

or

```julia
using Pkg; Pkg.add("https://github.com/islent/PhysicalParticles.jl")
```

To test the Package:
```julia
]test PhysicalParticles
```

## Usage

We reserve the normal point type for mathematical calculations,
its only difference from PhysicalVector is that normal points do not need definition of units.
To be brief we only show examples for PhysicalVector.

### Physical vectors

`Unitful` and `UnitfulAstro` are not necessary if you only use the inner defined vector types for simulation such as `Position` and `Velocity`.

1. Basic

    The default constructors always return a zero vector with the default unit:
    ```julia
    julia> PhysicalVector2D()
    PhysicalVector2D(0.0 m, 0.0 m)

    julia> PhysicalVector3D()
    PhysicalVector3D(0.0 m, 0.0 m, 0.0 m)
    ```

    Function `PhysicalVector` determines the dimension automatically, with units provided:
    ```julia
    julia> p = PhysicalVector(1.0,2.0,3.0,u"km")
    PhysicalVector3D(1.0 km, 2.0 km, 3.0 km)

    julia> v = PhysicalVector(1.0,2.0,3.0,u"km/s")
    PhysicalVector3D(1.0 km s^-1, 2.0 km s^-1, 3.0 km s^-1)
    ```

    Basic linear operators are overloaded and dimension-safe:
    ```julia
    julia> p+v
    ERROR: DimensionError: 1.0 km and 1.0 km s^-1 are not dimensionally compatible.

    julia> p + 2.0u"m"
    PhysicalVector3D(1002.0 m, 2002.0 m, 3002.0 m)

    julia> p + v * 10u"s"
    PhysicalVector3D(11.0 km, 22.0 km, 33.0 km)

    julia> 2p
    PhysicalVector3D(2.0 km, 4.0 km, 6.0 km)

    julia> p/4.0u"km/s"
    PhysicalVector3D(0.25 s, 0.5 s, 0.75 s)

    julia> norm(p)
    3.7416573867739413 km

    julia> dot(p, PhysicalVector(4.0,5.0,6.0,u"km"))
    32.0 km^2
    ```

2. Frequently used vector types

    For scientific simulations, we provide `Postition`, `PositionAstro`, `Velocity`, `VelocityAstro`, `Acceleration`, `AccelerationAstro` types.
    Vectors marked as `Astro` use `kpc` and `Gyr` as default unit, and those not markes use `m` and `s` as default.
    Here we display some basic operations:

    ```julia
    julia> p = PositionAstro(1.0,2.0,3.0)
    PhysicalVector3D(1.0 kpc, 2.0 kpc, 3.0 kpc)

    julia> v = VelocityAstro(1.0,1.0,1.0)
    PhysicalVector3D(1.0 kpc Gyr^-1, 1.0 kpc Gyr^-1, 1.0 kpc Gyr^-1)

    julia> a = AccelerationAstro(1.0,1.0,1.0)
    PhysicalVector3D(1.0 kpc Gyr^-2, 1.0 kpc Gyr^-2, 1.0 kpc Gyr^-2)

    julia> t = 0.1u"Gyr"
    0.1 Gyr

    julia> p2 = p + v*t + 0.5a*t^2
    PhysicalVector3D(1.105 kpc, 2.105 kpc, 3.105 kpc)
    ```

    Linear algebra is same with `PhysicalVector`

3. Manipulate arrays

    An example:
    ```julia
    julia> b = rand(3, 5)
    3×5 Array{Float64,2}:
    0.698576  0.873738  0.48515   0.216954   0.448364
    0.324348  0.474831  0.83997   0.0383103  0.317231
    0.296494  0.73641   0.494008  0.801979   0.686554

    julia> bv = pconvert(b, u"m")
    5-element Array{PhysicalVector3D,1}:
    PhysicalVector3D(0.698576245896922 m, 0.3243483035589918 m, 0.2964935693642543 m)
    PhysicalVector3D(0.8737376618851125 m, 0.47483139458960366 m, 0.7364098881041843 m)
    PhysicalVector3D(0.4851496800831139 m, 0.839970486136373 m, 0.4940081925624349 m)
    PhysicalVector3D(0.21695371153286058 m, 0.03831034395213351 m, 0.8019785383245348 m)
    PhysicalVector3D(0.4483641735559649 m, 0.31723087507351955 m, 0.686554138865415 m)

    julia> bn = npconvert(b)
    5-element Array{Point3D,1}:
    Point3D(0.698576245896922, 0.3243483035589918, 0.2964935693642543)
    Point3D(0.8737376618851125, 0.47483139458960366, 0.7364098881041843)
    Point3D(0.4851496800831139, 0.839970486136373, 0.4940081925624349)
    Point3D(0.21695371153286058, 0.03831034395213351, 0.8019785383245348)
    Point3D(0.4483641735559649, 0.31723087507351955, 0.686554138865415)

    julia> Center = center(bv)
    PhysicalVector3D(0.5453456867089865 m, 0.4391404150442533 m, 0.5492360538443946 m)

    julia> min_x(bv)
    0.21695371153286058 m
    ```

4. Physical Particles
    1. The default simulating particles are astronomical stars, whereas easy to be adapted in other simulation fields:
        ```julia
        julia> PhysicalParticle()
        PhysicalParticle3D(PhysicalVector3D(0.0 kpc, 0.0 kpc, 0.0 kpc), PhysicalVector3D(0.0 kpc Gyr^-1, 0.0 kpc Gyr^-1, 0.0 kpc Gyr^-1), PhysicalVector3D(0.0 kpc Gyr^-2, 0.0 kpc Gyr^-2, 0.0 kpc Gyr^-2), 0.0 M⊙, 0, star::ParticleType = 5)
        ```
        The struct contains the most basic information in numerical simulation:
        ```julia
        mutable struct PhysicalParticle3D <: AbstractParticle3D
            Pos::PhysicalVector3D
            Vel::PhysicalVector3D
            Acc::PhysicalVector3D
            Mass::Quantity
            ID::Int64
            Type::ParticleType
        end
        ```
        Particle types are defined in here:
        ```julia
        @enum ParticleType begin
            gas = 1
            halo = 2
            disk = 3
            bulge = 4
            star = 5
            blackhole = 6
        end
        ---------------------
        julia> ParticleType(1)
        gas::ParticleType = 1
        ```
    2. To cater for gas physics, there is a SPH (smoothed particle hydrodynamics) data type consistent with `Gadget2`
        ```julia
        mutable struct GasParticle2D <: AbstractParticle2D
            Pos::PhysicalVector2D
            Vel::PhysicalVector2D
            Acc::PhysicalVector2D
            Mass::Quantity
            ID::Int64

            Entropy::Quantity
            Density::Quantity
            Hsml::Quantity

            Left::Float64
            Right::Float64
            NumNgbFound::Int64

            RotVel::PhysicalVector2D
            DivVel::Quantity
            CurlVel::Quantity
            dHsmlRho::Float64

            Pressure::Quantity
            DtEntropy::Quantity
            MaxSignalVel::Quantity
        end
        ```
        You could even find the `Gadget2` header here:
        ```julia
        julia> Header_Gadget2()
        Header_Gadget2(Int32[0, 0, 0, 0, 0, 0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.0, 0.0, 0, 0, UInt32[0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000], 0, 1, 0.0, 0.3, 0.7, 0.71, 0, 0, UInt32[0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000], 0)
        ```
        And read `Gadget2` formatted snapshots by
        ```julia
        function read_gadget2(filename::String)
        ```
        which returns a tuple `(::Header_Gadget2, ::Array{PhysicalParticle,1}, ::Array{GasData,1})`
        Later on we would support more particle types, so the tuple would be expanded then.
        Now write a snapshot at will:
        ```julia
        function write_gadget2(filename::String, Header::Header_Gadget2, Particles::Array{PhysicalParticle,1}, SphData::Array{GasData,1})
        ```
5. Output an array of physical vectors or particles by simply calling function `write_ascii`
6. There is a physical constant struct containing the most useful constants in astrophysical simulations, supported by [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl):
    ```julia
    struct PhysicalConstant
        c::Constant # light speed
        G::Constant # Newtonian constant of gravitation
        h::Constant # Planck constant
        e::Constant # Elementary charge
        m_e::Constant # Electron mass
        m_n::Constant # Neutron mass
        m_p::Constant # Protron mass
        stefan_boltzmann::Constant # Stefan-Boltzmann constant
        H::Constant # Hubble constant

        ACC0::Constant # Modified gravitational acceleration constant
    end
    ```
    where H and ACC0 are user defined (not from CODATA2014):
    ```julia
    julia> Constants.H
    Hubble constant (H)
    Value                         = 74.03 km Mpc^-1 s^-1
    Standard uncertainty          = 1.42 km Mpc^-1 s^-1
    Relative standard uncertainty = 0.019
    Reference                     = Hubble Space Telescope 2019-03-18

    julia> Constants.ACC0
    Modified gravitational acceleration (ACC0)
    Value                         = 1.2e-8 cm s^-2
    Standard uncertainty          = (exact)
    Relative standard uncertainty = (exact)
    Reference                     = Milgrom 1983
    ```

4. Here are all of the exported types or functions:

    ```julia
    export
        AbstractPoint,
            AbstractPoint2D,
            AbstractPoint3D,
            Point, Point2D, Point3D,
            Position, Velocity, Acceleration,
            PositionAstro, VelocityAstro, AccelerationAstro,
            PhysicalVector, PhysicalVector2D, PhysicalVector3D,

            AbstractParticle, AbstractParticle2D, AbstractParticle3D,
            PhysicalParticle, PhysicalParticle2D, PhysicalParticle3D,

            ParticleType,
            Extent, Extent2D, Extent3D,
            PhysicalExtent, PhysicalExtent2D, PhysicalExtent3D,

            PhysicalConstant,

            getx, gety, getz,

            +,-,*,/,zero,length,iterate,

            # Linear Algebra
            norm, normalize, dot, cross,
            rotate, rotate_x, rotate_y, rotate_z,

            # Simulation Box
            mean,
            min_x, min_y, min_z,
            max_x, max_y, max_z,
            min_coord, max_coord,
            center_x, center_y, center_z,
            center,

            # Convert array to points
            pconvert, npconvert,

            # Peano-Hilbert algorithms
            peanokey, hilbertsort!, mssort!
    ```

### Physical Particles
