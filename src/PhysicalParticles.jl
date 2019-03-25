# PhysicalParticles.jl
# Author: islent (Runyu Meng, USTC, leoislent@gmail.com)
# Instructed by longqian95 (Qian Long, YNAO, longqian@ynao.ac.cn)
# References: GeometicalPredicates.jl and Gadget2

module PhysicalParticles

using Unitful, UnitfulAstro

import Unitful: Units
import Base: +,-,*,/,zero,length,iterate
import LinearAlgebra: norm, dot, cross

export
    AbstractPoint,
        AbstractPoint2D,
        AbstractPoint3D,
    Point, Point2D, Point3D,
    Pos, Vel, Acc,
    PosAstro, VelAstro, AccAstro,

    PhysicalVector, PhysicalVector2D, PhysicalVector3D,
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

############      Physical Vectors       ###########

struct PhysicalVector2D <: AbstractPoint3D
    x::Quantity
    y::Quantity
    PhysicalVector2D(x::Quantity,y::Quantity) = new(x, y)
    PhysicalVector2D(x::Real,y::Real, u::Units) = new(x*u, y*u)
end
PhysicalVector2D() = PhysicalVector2D(0.0u"m", 0.0u"m")

@inline getx(p::PhysicalVector2D) = p.x
@inline gety(p::PhysicalVector2D) = p.y

struct PhysicalVector3D <: AbstractPoint3D
    x::Quantity
    y::Quantity
    z::Quantity
    PhysicalVector3D(x::Quantity,y::Quantity,z::Quantity) = new(x, y, z)
    PhysicalVector3D(x::Real,y::Real,z::Real, u::Units) = new(x*u, y*u, z*u)
end
PhysicalVector3D() = PhysicalVector3D(0.0u"m", 0.0u"m", 0.0u"m")

@inline getx(p::PhysicalVector3D) = p.x
@inline gety(p::PhysicalVector3D) = p.y
@inline getz(p::PhysicalVector3D) = p.z

PhysicalVector(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
PhysicalVector(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)

############      Position       ###########

Pos(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
Pos(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
Pos(x::Real, y::Real) = PhysicalVector2D(x*u"m", y*u"m")
Pos(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"m", y*u"m", z*u"m")
PosAstro(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
PosAstro(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
PosAstro(x::Real, y::Real) = PhysicalVector2D(x*u"kpc", y*u"kpc")
PosAstro(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"kpc", y*u"kpc", z*u"kpc")

############      Position       ###########

Vel(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
Vel(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
Vel(x::Real, y::Real) = PhysicalVector2D(x*u"m/s", y*u"m/s")
Vel(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"m/s", y*u"m/s", z*u"m/s")
VelAstro(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
VelAstro(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
VelAstro(x::Real, y::Real) = PhysicalVector2D(x*u"kpc/Gyr", y*u"kpc/Gyr")
VelAstro(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"kpc/Gyr", y*u"kpc/Gyr", z*u"kpc/Gyr")

############      Position       ###########

Acc(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
Acc(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
Acc(x::Real, y::Real) = PhysicalVector2D(x*u"m/s^2", y*u"m/s^2")
Acc(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"m/s^2", y*u"m/s^2", z*u"m/s^2")
AccAstro(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
AccAstro(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
AccAstro(x::Real, y::Real) = PhysicalVector2D(x*u"kpc/Gyr^2", y*u"kpc/Gyr^2")
AccAstro(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"kpc/Gyr^2", y*u"kpc/Gyr^2", z*u"kpc/Gyr^2")

############      Basic mathematics       ###########

@inline norm(p::AbstractPoint2D) = sqrt(p.x^2 + p.y^2)
@inline norm(p::AbstractPoint3D) = sqrt(p.x^2 + p.y^2 + p.z^2)
@inline length(p::AbstractPoint) = 1
@inline iterate(p::AbstractPoint) = (p,nothing)
@inline iterate(p::AbstractPoint,st) = nothing
@inline real(p::AbstractPoint) = p

# 2D
@inline +(p1::PhysicalVector2D, p2::PhysicalVector2D) = PhysicalVector2D(p1.x+p2.x, p1.y+p2.y)
@inline -(p1::PhysicalVector2D, p2::PhysicalVector2D) = PhysicalVector2D(p1.x-p2.x, p1.y-p2.y)
@inline *(p1::PhysicalVector2D, p2::PhysicalVector2D) = PhysicalVector2D(p1.x*p2.x + p1.y*p2.y)
@inline +(p1::PhysicalVector2D, a::Quantity) = PhysicalVector2D(p1.x+a, p1.y+a)
@inline -(p1::PhysicalVector2D, a::Quantity) = PhysicalVector2D(p1.x-a, p1.y-a)
@inline *(p1::PhysicalVector2D, a::Quantity) = PhysicalVector2D(p1.x*a, p1.y*a)
@inline /(p1::PhysicalVector2D, a::Quantity) = PhysicalVector2D(p1.x/a, p1.y/a)
@inline +(a::Quantity, p1::PhysicalVector2D) = PhysicalVector2D(p1.x+a, p1.y+a)
@inline -(a::Quantity, p1::PhysicalVector2D) = PhysicalVector2D(p1.x-a, p1.y-a)
@inline *(a::Quantity, p1::PhysicalVector2D) = PhysicalVector2D(p1.x*a, p1.y*a)
@inline /(a::Quantity, p1::PhysicalVector2D) = PhysicalVector2D(p1.x/a, p1.y/a)
@inline dot(p1::PhysicalVector2D, p2::PhysicalVector2D) = PhysicalVector2D(p1.x*p2.x + p1.y*p2.y)

# 3D
@inline +(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z)
@inline -(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
@inline *(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z)
@inline +(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x+a, p1.y+a, p1.z+a)
@inline -(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x-a, p1.y-a, p1.z-a)
@inline *(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x*a, p1.y*a, p1.z*a)
@inline /(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x/a, p1.y/a, p1.z/a)
@inline +(a::Quantity, p1::PhysicalVector3D) = PhysicalVector3D(p1.x+a, p1.y+a, p1.z+a)
@inline -(a::Quantity, p1::PhysicalVector3D) = PhysicalVector3D(p1.x-a, p1.y-a, p1.z-a)
@inline *(a::Quantity, p1::PhysicalVector3D) = PhysicalVector3D(p1.x*a, p1.y*a, p1.z*a)
@inline /(a::Quantity, p1::PhysicalVector3D) = PhysicalVector3D(p1.x/a, p1.y/a, p1.z/a)
@inline dot(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z)

############      Linear Algebra       ###########



############      Physical Particles       ###########


end  # module PhysicalParticles
