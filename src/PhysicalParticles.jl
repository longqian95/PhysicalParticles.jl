# PhysicalParticles.jl
# Author: islent (Runyu Meng, USTC, leoislent@gmail.com)
# Instructed by longqian95 (Qian Long, YNAO, longqian@ynao.ac.cn)
# References: GeometicalPredicates.jl and Gadget2

module PhysicalParticles

using Unitful, UnitfulAstro

import Base: +,-,*,/,zero,length,iterate
import LinearAlgebra: norm, dot, cross

export
    AbstractPoint,
        AbstractPoint2D,
        AbstractPoint3D,
    Point, Point2D, Point3D,
    Pos, Vel, Acc,
    PosAstro, VelAstro, AccAstro,
    Pos2D, Vel2D, Acc2D,
    Pos3D, Vel3D, Acc3D,
    PosAstro2D, VelAstro2D, AccAstro2D,
    PosAstro3D, VelAstro3D, AccAstro3D,

    PhysicalParticle, PhysicalParticleAstro,

    getx, gety, getz,


    norm

    #peanokey, hilbertsort!


abstract type AbstractPoint end
abstract type AbstractPoint2D <: AbstractPoint end
abstract type AbstractPoint3D <: AbstractPoint end

# standard 2D point
struct Point2D <: AbstractPoint2D
    x::Float64
    y::Float64
    Point2D(x::Float64,y::Float64) = new(x, y)
end
Point2D() = Point2D(0.0, 0.0)

@inline getx(p::Point2D) = p.x
@inline gety(p::Point2D) = p.y

# standard 3D point
struct Point3D <: AbstractPoint3D
    x::Float64
    y::Float64
    z::Float64
    Point3D(x::Float64,y::Float64,z::Float64) = new(x, y, z)
end
Point3D() = Point3D(0.0, 0.0, 0.0)

@inline getx(p::Point3D) = p.x
@inline gety(p::Point3D) = p.y
@inline getz(p::Point3D) = p.z

Point(x::Real, y::Real) = Point2D(x, y)
Point(x::Real, y::Real, z::Real) = Point3D(x, y, z)

############      Position       ###########
# standard 2D Position
struct Pos2D <: AbstractPoint2D
    x::Quantity
    y::Quantity
    Pos2D(x::Quantity,y::Quantity) = new(x, y)
    Pos2D(x::Real,y::Real) = new(x*u"m", y*u"m")
end
Pos2D() = Pos2D(0.0u"m", 0.0u"m")

@inline getx(p::Pos2D) = p.x
@inline gety(p::Pos2D) = p.y

# standard 3D Position
struct Pos3D <: AbstractPoint3D
    x::Quantity
    y::Quantity
    z::Quantity
    Pos3D(x::Quantity,y::Quantity,z::Quantity) = new(x, y, z)
    Pos3D(x::Real,y::Real,z::Real) = new(x*u"m", y*u"m", z*u"m")
end
Pos3D() = Pos3D(0.0u"m", 0.0u"m", 0.0u"m")

@inline getx(p::Pos3D) = p.x
@inline gety(p::Pos3D) = p.y
@inline getz(p::Pos3D) = p.z

Pos(x::Real, y::Real) = Pos2D(x*u"m", y*u"m")
Pos(x::Real, y::Real, z::Real) = Pos3D(x*u"m", y*u"m", z*u"m")

# Astro 2D Position
struct PosAstro2D <: AbstractPoint2D
    x::Quantity
    y::Quantity
    PosAstro2D(x::Quantity,y::Quantity) = new(x, y)
    PosAstro2D(x::Real,y::Real) = new(x*u"kpc", y*u"kpc")
end
PosAstro2D() = PosAstro2D(0.0u"kpc", 0.0u"kpc")

@inline getx(p::PosAstro2D) = p.x
@inline gety(p::PosAstro2D) = p.y

# Astro 3D Position
struct PosAstro3D <: AbstractPoint3D
    x::Quantity
    y::Quantity
    z::Quantity
    PosAstro3D(x::Quantity,y::Quantity,z::Quantity) = new(x, y, z)
    PosAstro3D(x::Real,y::Real,z::Real) = new(x*u"kpc", y*u"kpc", z*u"kpc")
end
PosAstro3D() = PosAstro3D(0.0u"kpc", 0.0u"kpc", 0.0u"kpc")

@inline getx(p::PosAstro3D) = p.x
@inline gety(p::PosAstro3D) = p.y
@inline getz(p::PosAstro3D) = p.z

PosAstro(x::Real, y::Real) = PosAstro2D(x*u"kpc", y*u"kpc")
PosAstro(x::Real, y::Real, z::Real) = PosAstro3D(x*u"kpc", y*u"kpc", z*u"kpc")

############      Velocity       ###########
# standard 2D Velocity
struct Vel2D <: AbstractPoint2D
    x::Quantity
    y::Quantity
    Vel2D(x::Quantity,y::Quantity) = new(x, y)
    Vel2D(x::Real,y::Real) = new(x*u"m/s", y*u"m/s")
end
Vel2D() = Vel2D(0.0u"m/s", 0.0u"m/s")

@inline getx(p::Vel2D) = p.x
@inline gety(p::Vel2D) = p.y

# standard 3D Velocity
struct Vel3D <: AbstractPoint3D
    x::Quantity
    y::Quantity
    z::Quantity
    Vel3D(x::Quantity,y::Quantity,z::Quantity) = new(x, y, z)
    Vel3D(x::Real,y::Real,z::Real) = new(x*u"m/s", y*u"m/s", z*u"m/s")
end
Vel3D() = Vel3D(0.0u"m/s", 0.0u"m/s", 0.0u"m/s")

@inline getx(p::Vel3D) = p.x
@inline gety(p::Vel3D) = p.y
@inline getz(p::Vel3D) = p.z

Vel(x::Real, y::Real) = Vel2D(x*u"m/s", y*u"m/s")
Vel(x::Real, y::Real, z::Real) = Vel3D(x*u"m/s", y*u"m/s", z*u"m/s")

# Astro 2D Velocity
struct VelAstro2D <: AbstractPoint2D
    x::Quantity
    y::Quantity
    VelAstro2D(x::Quantity,y::Quantity) = new(x, y)
    VelAstro2D(x::Real,y::Real) = new(x*u"kpc/Gyr", y*u"kpc/Gyr")
end
VelAstro2D() = VelAstro2D(0.0u"kpc/Gyr", 0.0u"kpc/Gyr")

@inline getx(p::VelAstro2D) = p.x
@inline gety(p::VelAstro2D) = p.y

# Astro 3D Velocity
struct VelAstro3D <: AbstractPoint3D
    x::Quantity
    y::Quantity
    z::Quantity
    VelAstro3D(x::Quantity,y::Quantity,z::Quantity) = new(x, y, z)
    VelAstro3D(x::Real,y::Real,z::Real) = new(x*u"kpc/Gyr", y*u"kpc/Gyr", z*u"kpc/Gyr")
end
VelAstro3D() = VelAstro3D(0.0u"kpc/Gyr", 0.0u"kpc/Gyr", 0.0u"kpc/Gyr")

@inline getx(p::VelAstro3D) = p.x
@inline gety(p::VelAstro3D) = p.y
@inline getz(p::VelAstro3D) = p.z

VelAstro(x::Real, y::Real) = VelAstro2D(x*u"kpc/Gyr", y*u"kpc/Gyr")
VelAstro(x::Real, y::Real, z::Real) = VelAstro3D(x*u"kpc/Gyr", y*u"kpc/Gyr", z*u"kpc/Gyr")

############      Acceleration       ###########
# standard 2D Acceleration
struct Acc2D <: AbstractPoint2D
    x::Quantity
    y::Quantity
    Acc2D(x::Quantity,y::Quantity) = new(x, y)
    Acc2D(x::Real,y::Real) = new(x*u"m/s^2", y*u"m/s^2")
end
Acc2D() = Acc2D(0.0u"m/s^2", 0.0u"m/s^2")

@inline getx(p::Acc2D) = p.x
@inline gety(p::Acc2D) = p.y

# standard 3D Acceleration
struct Acc3D <: AbstractPoint3D
    x::Quantity
    y::Quantity
    z::Quantity
    Acc3D(x::Quantity,y::Quantity,z::Quantity) = new(x, y, z)
    Acc3D(x::Real,y::Real,z::Real) = new(x*u"m/s^2", y*u"m/s^2", z*u"m/s^2")
end
Acc3D() = Acc3D(0.0u"m/s^2", 0.0u"m/s^2", 0.0u"m/s^2")

@inline getx(p::Acc3D) = p.x
@inline gety(p::Acc3D) = p.y
@inline getz(p::Acc3D) = p.z

Acc(x::Real, y::Real) = Acc2D(x*u"m/s^2", y*u"m/s^2")
Acc(x::Real, y::Real, z::Real) = Acc3D(x*u"m/s^2", y*u"m/s^2", z*u"m/s^2")

# Astro 2D Acceleration
struct AccAstro2D <: AbstractPoint2D
    x::Quantity
    y::Quantity
    AccAstro2D(x::Quantity,y::Quantity) = new(x, y)
    AccAstro2D(x::Real,y::Real) = new(x*u"kpc/Gyr^2", y*u"kpc/Gyr^2")
end
AccAstro2D() = AccAstro2D(0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2")

@inline getx(p::AccAstro2D) = p.x
@inline gety(p::AccAstro2D) = p.y

# Astro 3D Acceleration
struct AccAstro3D <: AbstractPoint3D
    x::Quantity
    y::Quantity
    z::Quantity
    AccAstro3D(x::Quantity,y::Quantity,z::Quantity) = new(x, y, z)
    AccAstro3D(x::Real,y::Real,z::Real) = new(x*u"kpc/Gyr^2", y*u"kpc/Gyr^2", z*u"kpc/Gyr^2")
end
AccAstro3D() = AccAstro3D(0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2")

@inline getx(p::AccAstro3D) = p.x
@inline gety(p::AccAstro3D) = p.y
@inline getz(p::AccAstro3D) = p.z

AccAstro(x::Real, y::Real) = AccAstro2D(x*u"kpc/Gyr^2", y*u"kpc/Gyr^2")
AccAstro(x::Real, y::Real, z::Real) = AccAstro3D(x*u"kpc/Gyr^2", y*u"kpc/Gyr^2", z*u"kpc/Gyr^2")

############      Basic mathematics       ###########

@inline norm(p::AbstractPoint2D) = sqrt(p.x^2 + p.y^2)
@inline norm(p::AbstractPoint3D) = sqrt(p.x^2 + p.y^2 + p.z^2)
@inline length(p::AbstractPoint) = 1
@inline iterate(p::AbstractPoint) = (p,nothing)
@inline iterate(p::AbstractPoint,st) = nothing
@inline real(p::AbstractPoint) = p

@inline +(p1::Point2D, p2::Point2D) = Point2D(p1.x+p2.x, p1.y+p2.y)
@inline -(p1::Point2D, p2::Point2D) = Point2D(p1.x-p2.x, p1.y-p2.y)
@inline *(p1::Point2D, p2::Point2D) = Point2D(p1.x*p2.x + p1.y*p2.y)

@inline +(p1::Point3D, p2::Point3D) = Point3D(p1.x+p2.x, p1.y+p2.y, p1.y+p2.y)
@inline -(p1::Point3D, p2::Point3D) = Point3D(p1.x-p2.x, p1.y-p2.y, p1.y-p2.y)
@inline *(p1::Point3D, p2::Point3D) = Point3D(p1.x*p2.x + p1.y*p2.y + p1.y*p2.y)

@inline +(p1::Pos2D, p2::Pos2D) = Pos2D(p1.x+p2.x, p1.y+p2.y)
@inline -(p1::Pos2D, p2::Pos2D) = Pos2D(p1.x-p2.x, p1.y-p2.y)
@inline *(p1::Pos2D, p2::Pos2D) = Pos2D(p1.x*p2.x + p1.y*p2.y)

@inline +(p1::Pos3D, p2::Pos3D) = Pos3D(p1.x+p2.x, p1.y+p2.y, p1.y+p2.y)
@inline -(p1::Pos3D, p2::Pos3D) = Pos3D(p1.x-p2.x, p1.y-p2.y, p1.y-p2.y)
@inline *(p1::Pos3D, p2::Pos3D) = Pos3D(p1.x*p2.x + p1.y*p2.y + p1.y*p2.y)

@inline +(p1::Vel2D, p2::Vel2D) = Vel2D(p1.x+p2.x, p1.y+p2.y)
@inline -(p1::Vel2D, p2::Vel2D) = Vel2D(p1.x-p2.x, p1.y-p2.y)
@inline *(p1::Vel2D, p2::Vel2D) = Vel2D(p1.x*p2.x + p1.y*p2.y)

@inline +(p1::Vel3D, p2::Vel3D) = Vel3D(p1.x+p2.x, p1.y+p2.y, p1.y+p2.y)
@inline -(p1::Vel3D, p2::Vel3D) = Vel3D(p1.x-p2.x, p1.y-p2.y, p1.y-p2.y)
@inline *(p1::Vel3D, p2::Vel3D) = Vel3D(p1.x*p2.x + p1.y*p2.y + p1.y*p2.y)

@inline +(p1::Acc2D, p2::Acc2D) = Acc2D(p1.x+p2.x, p1.y+p2.y)
@inline -(p1::Acc2D, p2::Acc2D) = Acc2D(p1.x-p2.x, p1.y-p2.y)
@inline *(p1::Acc2D, p2::Acc2D) = Acc2D(p1.x*p2.x + p1.y*p2.y)

@inline +(p1::Acc3D, p2::Acc3D) = Acc3D(p1.x+p2.x, p1.y+p2.y, p1.y+p2.y)
@inline -(p1::Acc3D, p2::Acc3D) = Acc3D(p1.x-p2.x, p1.y-p2.y, p1.y-p2.y)
@inline *(p1::Acc3D, p2::Acc3D) = Acc3D(p1.x*p2.x + p1.y*p2.y + p1.y*p2.y)

@inline +(p1::PosAstro2D, p2::PosAstro2D) = PosAstro2D(p1.x+p2.x, p1.y+p2.y)
@inline -(p1::PosAstro2D, p2::PosAstro2D) = PosAstro2D(p1.x-p2.x, p1.y-p2.y)
@inline *(p1::PosAstro2D, p2::PosAstro2D) = PosAstro2D(p1.x*p2.x + p1.y*p2.y)

@inline +(p1::PosAstro3D, p2::PosAstro3D) = PosAstro3D(p1.x+p2.x, p1.y+p2.y, p1.y+p2.y)
@inline -(p1::PosAstro3D, p2::PosAstro3D) = PosAstro3D(p1.x-p2.x, p1.y-p2.y, p1.y-p2.y)
@inline *(p1::PosAstro3D, p2::PosAstro3D) = PosAstro3D(p1.x*p2.x + p1.y*p2.y + p1.y*p2.y)

@inline +(p1::VelAstro2D, p2::VelAstro2D) = VelAstro2D(p1.x+p2.x, p1.y+p2.y)
@inline -(p1::VelAstro2D, p2::VelAstro2D) = VelAstro2D(p1.x-p2.x, p1.y-p2.y)
@inline *(p1::VelAstro2D, p2::VelAstro2D) = VelAstro2D(p1.x*p2.x + p1.y*p2.y)

@inline +(p1::VelAstro3D, p2::VelAstro3D) = VelAstro3D(p1.x+p2.x, p1.y+p2.y, p1.y+p2.y)
@inline -(p1::VelAstro3D, p2::VelAstro3D) = VelAstro3D(p1.x-p2.x, p1.y-p2.y, p1.y-p2.y)
@inline *(p1::VelAstro3D, p2::VelAstro3D) = VelAstro3D(p1.x*p2.x + p1.y*p2.y + p1.y*p2.y)

@inline +(p1::AccAstro2D, p2::AccAstro2D) = AccAstro2D(p1.x+p2.x, p1.y+p2.y)
@inline -(p1::AccAstro2D, p2::AccAstro2D) = AccAstro2D(p1.x-p2.x, p1.y-p2.y)
@inline *(p1::AccAstro2D, p2::AccAstro2D) = AccAstro2D(p1.x*p2.x + p1.y*p2.y)

@inline +(p1::AccAstro3D, p2::AccAstro3D) = AccAstro3D(p1.x+p2.x, p1.y+p2.y, p1.y+p2.y)
@inline -(p1::AccAstro3D, p2::AccAstro3D) = AccAstro3D(p1.x-p2.x, p1.y-p2.y, p1.y-p2.y)
@inline *(p1::AccAstro3D, p2::AccAstro3D) = AccAstro3D(p1.x*p2.x + p1.y*p2.y + p1.y*p2.y)

############      Linear Algebra       ###########

############      Physical Particles       ###########


end  # module PhysicalParticles
