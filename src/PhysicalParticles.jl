# PhysicalParticles.jl
# Author: islent (Runyu Meng, USTC, leoislent@gmail.com)
# Instructed by longqian95 (Qian Long, YNAO, longqian@ynao.ac.cn)
# References: GeometicalPredicates.jl and Gadget2

module PhysicalParticles

using Unitful, UnitfulAstro

import Unitful: Units
import Base: +,-,*,/,zero,length,iterate,
            rand
import LinearAlgebra: norm, normalize, normalize!, dot, cross

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

    +,-,*,/,zero,length,iterate,


    norm, normalize, normalize!, dot, cross,

    mean,
    pconvert
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

struct PhysicalVector2D <: AbstractPoint2D
    x::Quantity
    y::Quantity
    PhysicalVector2D(x::Quantity,y::Quantity) = new(x, y)
    PhysicalVector2D(x::Real,y::Real, u::Units) = new(x*u, y*u)
end
PhysicalVector2D() = PhysicalVector2D(0.0u"m", 0.0u"m")
PhysicalVector2D(u::Units) = PhysicalVector2D(0.0u, 0.0u)

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
PhysicalVector3D(u::Units) = PhysicalVector3D(0.0u, 0.0u, 0.0u)

@inline getx(p::PhysicalVector3D) = p.x
@inline gety(p::PhysicalVector3D) = p.y
@inline getz(p::PhysicalVector3D) = p.z

PhysicalVector(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
PhysicalVector(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)

############      Position       ###########

Pos() = PhysicalVector3D(0.0u"m", 0.0u"m", 0.0u"m")
Pos(u::Units) = PhysicalVector3D(0.0u, 0.0u, 0.0u)
Pos(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
Pos(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
Pos(x::Real, y::Real) = PhysicalVector2D(x*u"m", y*u"m")
Pos(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"m", y*u"m", z*u"m")

PosAstro() = PhysicalVector3D(0.0u"kpc", 0.0u"kpc", 0.0u"kpc")
PosAstro(u::Units) = PhysicalVector3D(0.0u, 0.0u, 0.0u)
PosAstro(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
PosAstro(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
PosAstro(x::Real, y::Real) = PhysicalVector2D(x*u"kpc", y*u"kpc")
PosAstro(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"kpc", y*u"kpc", z*u"kpc")

############      Position       ###########

Vel() = PhysicalVector3D(0.0u"m/s", 0.0u"m/s", 0.0u"m/s")
Vel(u::Units) = PhysicalVector3D(0.0u, 0.0u, 0.0u)
Vel(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
Vel(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
Vel(x::Real, y::Real) = PhysicalVector2D(x*u"m/s", y*u"m/s")
Vel(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"m/s", y*u"m/s", z*u"m/s")

VelAstro() = PhysicalVector3D(0.0u"kpc/Gyr", 0.0u"kpc/Gyr", 0.0u"kpc/Gyr")
VelAstro(u::Units) = PhysicalVector3D(0.0u, 0.0u, 0.0u)
VelAstro(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
VelAstro(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
VelAstro(x::Real, y::Real) = PhysicalVector2D(x*u"kpc/Gyr", y*u"kpc/Gyr")
VelAstro(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"kpc/Gyr", y*u"kpc/Gyr", z*u"kpc/Gyr")

############      Position       ###########

Acc() = PhysicalVector3D(0.0u"m/s^2", 0.0u"m/s^2", 0.0u"m/s^2")
Acc(u::Units) = PhysicalVector3D(0.0u, 0.0u, 0.0u)
Acc(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
Acc(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
Acc(x::Real, y::Real) = PhysicalVector2D(x*u"m/s^2", y*u"m/s^2")
Acc(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"m/s^2", y*u"m/s^2", z*u"m/s^2")

AccAstro() = PhysicalVector3D(0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2")
AccAstro(u::Units) = PhysicalVector3D(0.0u, 0.0u, 0.0u)
AccAstro(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
AccAstro(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
AccAstro(x::Real, y::Real) = PhysicalVector2D(x*u"kpc/Gyr^2", y*u"kpc/Gyr^2")
AccAstro(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"kpc/Gyr^2", y*u"kpc/Gyr^2", z*u"kpc/Gyr^2")

############      Basic mathematics       ###########
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
@inline *(p1::PhysicalVector2D, a::Real) = PhysicalVector2D(p1.x*a, p1.y*a)
@inline /(p1::PhysicalVector2D, a::Real) = PhysicalVector2D(p1.x/a, p1.y/a)
@inline +(a::Quantity, p1::PhysicalVector2D) = PhysicalVector2D(p1.x+a, p1.y+a)
@inline -(a::Quantity, p1::PhysicalVector2D) = PhysicalVector2D(p1.x-a, p1.y-a)
@inline *(a::Quantity, p1::PhysicalVector2D) = PhysicalVector2D(p1.x*a, p1.y*a)
@inline /(a::Quantity, p1::PhysicalVector2D) = PhysicalVector2D(p1.x/a, p1.y/a)
@inline *(a::Real, p1::PhysicalVector2D) = PhysicalVector2D(p1.x*a, p1.y*a)
@inline /(a::Real, p1::PhysicalVector2D) = PhysicalVector2D(p1.x/a, p1.y/a)
@inline norm(p::AbstractPoint2D) = sqrt(p.x^2 + p.y^2)
@inline dot(p1::PhysicalVector2D, p2::PhysicalVector2D) = PhysicalVector2D(p1.x*p2.x + p1.y*p2.y)
@inline zero(p::PhysicalVector2D) = PhysicalVector2D(p.x*0.0, p.y*0.0)

# 3D
@inline +(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z)
@inline -(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
@inline *(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z)
@inline *(p1::PhysicalVector3D, a::Real) = PhysicalVector3D(p1.x*a, p1.y*a, p1.z*a)
@inline +(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x+a, p1.y+a, p1.z+a)
@inline -(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x-a, p1.y-a, p1.z-a)
@inline *(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x*a, p1.y*a, p1.z*a)
@inline /(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x/a, p1.y/a, p1.z/a)
@inline *(p1::PhysicalVector3D, a::Real) = PhysicalVector3D(p1.x*a, p1.y*a, p1.z*a)
@inline /(p1::PhysicalVector3D, a::Real) = PhysicalVector3D(p1.x/a, p1.y/a, p1.z/a)
@inline +(a::Quantity, p1::PhysicalVector3D) = PhysicalVector3D(p1.x+a, p1.y+a, p1.z+a)
@inline -(a::Quantity, p1::PhysicalVector3D) = PhysicalVector3D(p1.x-a, p1.y-a, p1.z-a)
@inline *(a::Quantity, p1::PhysicalVector3D) = PhysicalVector3D(p1.x*a, p1.y*a, p1.z*a)
@inline /(a::Quantity, p1::PhysicalVector3D) = PhysicalVector3D(p1.x/a, p1.y/a, p1.z/a)
@inline *(a::Real, p1::PhysicalVector3D) = PhysicalVector3D(p1.x*a, p1.y*a, p1.z*a)
@inline /(a::Real, p1::PhysicalVector3D) = PhysicalVector3D(p1.x/a, p1.y/a, p1.z/a)
@inline norm(p::AbstractPoint3D) = sqrt(p.x^2 + p.y^2 + p.z^2)
@inline dot(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z)
@inline zero(p::PhysicalVector3D) = PhysicalVector3D(p.x*0.0, p.y*0.0, p.z*0.0)
@inline cross(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.y*p2.z-p1.z*p2.y, p1.z*p2.x-p1.x*p2.z, p1.x*p2.y-p1.y*p2.x)

############      Linear Algebra       ###########

@inline norm(p::PhysicalVector2D) = sqrt(p.x*p.x + p.y*p.y + p.z*p.z)
@inline normalize(p::PhysicalVector2D) = (n = norm(p); return PhysicalVector(p.x/n, p.y/n, p.z/n))
@inline normalize!(p::PhysicalVector2D) = (n = norm(p); p.x /= n; p.y /= n; p.z /= n)

function mean(a::Array{PhysicalVector3D})
    len = length(a)
    p = PhysicalVector3D()
    for i in 1:len
        p += a[i]
    end
    return p/len
end

function mean(a::Array{PhysicalVector2D})
    len = length(a)
    p = PhysicalVector3D()
    for i in 1:len
        p += a[i]
    end
    return p/len
end

function pconvert(a::Array{Float64,1})
    if length(a) == 3
        return PhysicalVector3D(a[1], a[2], a[3], u"m")
    elseif length(a) == 2
        return PhysicalVector(a[1], a[2], u"m")
    else
        error("Not supported dimension!")
    end
end

function pconvert(a::Array{Float64,1}, u)
    if length(a) == 3
        return PhysicalVector3D(a[1], a[2], a[3], u)
    elseif length(a) == 2
        return PhysicalVector(a[1], a[2], u)
    else
        error("Not supported dimension!")
    end
end

function pconvert(a::Array{Float64,2}, u)
    row, col = size(a)
    if row == 3
        p = rand(PhysicalVector3D,0)
        for i in 1:col
            @inbounds append!(p, PhysicalVector3D(a[1,i], a[2,i], a[3,i], u))
        end
        return p
    elseif row == 2
        p = rand(PhysicalVector2D,0)
        for i in 1:col
            @inbounds append!(p, PhysicalVector3D(a[1,i], a[2,i], u))
        end
        return p
    else
        error("Not supported dimension!")
    end
end

############      Physical Particles       ###########


end  # module PhysicalParticles
