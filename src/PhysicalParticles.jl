# PhysicalParticles.jl
# Author: islent (Runyu Meng, USTC, leoislent@gmail.com)
# Instructed by longqian95 (Qian Long, YNAO, longqian@ynao.ac.cn)
# References: GeometicalPredicates.jl and Gadget2

module PhysicalParticles

using Unitful, UnitfulAstro

import Unitful: Units
import Base: +,-,*,/,zero,length,iterate,
            rand
import LinearAlgebra: norm, normalize, dot, cross
import PhysicalConstants: CODATA2014, Constant, @constant

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

    GasParticle, GasParticle2D, GasParticle3D,
    GasData, GasData2D, GasData3D,

    ParticleType,
    Extent, Extent2D, Extent3D,
    PhysicalExtent, PhysicalExtent2D, PhysicalExtent3D,

    Constants,

    # Serve for ISLENT project
    PhysicalConstant,
    Header_Gadget2,
    TreeNode, PhysicalTreeNode,

    getx, gety, getz,

    +,-,*,/,zero,length,iterate,

    # Linear Algebra
    norm, normalize, dot, cross,
    rotate, rotate_x, rotate_y, rotate_z,
    distance,

    # Simulation Box
    mean,
    min_x, min_y, min_z,
    max_x, max_y, max_z,
    min_coord, max_coord,
    center_x, center_y, center_z,
    center, mass_center,

    # Convert array to points
    pconvert, npconvert,

    # Peano-Hilbert algorithms
    peanokey, hilbertsort!, mssort!,

    # file I/O
    write_ascii, read_ascii,
    write_gadget2, read_gadget2,

    # Numerics
    kdtree_setup, kdtree_k_search, kdtree_radius_search


abstract type AbstractPoint end
abstract type AbstractPoint2D <: AbstractPoint end
abstract type AbstractPoint3D <: AbstractPoint end

# We support normal support for mathematical applications
include("./NormalPoints.jl")

############      Physical Vectors       ###########
"""
    struct PhysicalVector2D <: AbstractPoint2D

        Fields
        ≡≡≡≡≡≡≡≡

        x :: Quantity
        y :: Quantity
"""
struct PhysicalVector2D <: AbstractPoint2D
    x::Quantity
    y::Quantity
    PhysicalVector2D(x::Quantity,y::Quantity) = new(x, y)
    PhysicalVector2D(x::Real,y::Real, u::Units=u"m") = new(x*u, y*u)
end

"""
    Constructors of struct PhysicalVector2D <: AbstractPoint2D
        Default unit: "m"

        PhysicalVector2D() -------------- Returns a zero vector
        PhysicalVector2D(u::Units) ------ Returns a zero vector
        PhysicalVector2D(x::Quantity,y::Quantity)
        PhysicalVector2D(x::Real,y::Real, u::Units=u"m")

        Examples
        ≡≡≡≡≡≡≡≡≡≡

        julia> PhysicalVector2D()
        PhysicalVector2D(0.0 m, 0.0 m)

        julia> a = PhysicalVector2D(1.0, 2.0)
        PhysicalVector2D(1.0 m, 2.0 m)

        julia> b = PhysicalVector2D(3,4,u"km")
        PhysicalVector2D(3 km, 4 km)

        julia> a+b
        PhysicalVector2D(3001.0 m, 4002.0 m)

        julia> norm(a)
        2.23606797749979 m
"""
PhysicalVector2D() = PhysicalVector2D(0.0u"m", 0.0u"m")
PhysicalVector2D(u::Units) = PhysicalVector2D(0.0u, 0.0u)

@inline getx(p::PhysicalVector2D) = p.x
@inline gety(p::PhysicalVector2D) = p.y

"""
    struct PhysicalVector3D <: AbstractPoint3D

        Fields
        ≡≡≡≡≡≡≡≡

        x :: Quantity
        y :: Quantity
        z :: Quantity
"""
struct PhysicalVector3D <: AbstractPoint3D
    x::Quantity
    y::Quantity
    z::Quantity
    PhysicalVector3D(x::Quantity,y::Quantity,z::Quantity) = new(x, y, z)
    PhysicalVector3D(x::Real,y::Real,z::Real, u::Units=u"m") = new(x*u, y*u, z*u)
end

"""
    Constructors of struct PhysicalVector3D <: AbstractPoint3D
        Default unit: "m"

        PhysicalVector3D() -------------- Returns a zero vector
        PhysicalVector3D(u::Units) ------ Returns a zero vector
        PhysicalVector3D(x::Quantity,y::Quantity,z::Quantity)
        PhysicalVector3D(x::Real,y::Real,z::Real, u::Units=u"m")

        Examples
        ≡≡≡≡≡≡≡≡≡≡

        julia> PhysicalVector3D()
        PhysicalVector3D(0.0 m, 0.0 m, 0.0 m)

        julia> a = PhysicalVector3D(1.0,2.0,3.0)
        PhysicalVector3D(1.0 m, 2.0 m, 3.0 m)

        julia> b = PhysicalVector3D(4.0,5.0,6.0, u"kpc")
        PhysicalVector3D(4.0 kpc, 5.0 kpc, 6.0 kpc)

        julia> a+b
        PhysicalVector3D(1.2342710325965468e20 m, 1.5428387907456834e20 m, 1.8514065488948203e20 m)
"""
PhysicalVector3D() = PhysicalVector3D(0.0u"m", 0.0u"m", 0.0u"m")
PhysicalVector3D(u::Units) = PhysicalVector3D(0.0u, 0.0u, 0.0u)

@inline getx(p::PhysicalVector3D) = p.x
@inline gety(p::PhysicalVector3D) = p.y
@inline getz(p::PhysicalVector3D) = p.z

"""
    Returns physical vectors

    PhysicalVector(x::Quantity, y::Quantity) = PhysicalVector2D(x,y)
    PhysicalVector(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
    PhysicalVector(x::Quantity, y::Quantity, z::Quantity) = PhysicalVector3D(x,y,z)
    PhysicalVector(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)
"""
PhysicalVector(x::Quantity, y::Quantity) = PhysicalVector2D(x,y)
PhysicalVector(x::Real, y::Real, u::Units) = PhysicalVector2D(x*u, y*u)
PhysicalVector(x::Quantity, y::Quantity, z::Quantity) = PhysicalVector3D(x,y,z)
PhysicalVector(x::Real, y::Real, z::Real, u::Units) = PhysicalVector3D(x*u, y*u, z*u)

############      Position       ###########
Position() = PhysicalVector3D(0.0u"m", 0.0u"m", 0.0u"m")
Position(u::Units=u"m") = PhysicalVector3D(0.0u, 0.0u, 0.0u)
Position(x::Real, y::Real, u::Units=u"m") = PhysicalVector2D(x*u, y*u)
Position(x::Real, y::Real, z::Real, u::Units=u"m") = PhysicalVector3D(x*u, y*u, z*u)
Position(x::Real, y::Real) = PhysicalVector2D(x*u"m", y*u"m")
Position(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"m", y*u"m", z*u"m")

PositionAstro() = PhysicalVector3D(0.0u"kpc", 0.0u"kpc", 0.0u"kpc")
PositionAstro(u::Units=u"kpc") = PhysicalVector3D(0.0u, 0.0u, 0.0u)
PositionAstro(x::Real, y::Real, u::Units=u"kpc") = PhysicalVector2D(x*u, y*u)
PositionAstro(x::Real, y::Real, z::Real, u::Units=u"kpc") = PhysicalVector3D(x*u, y*u, z*u)
PositionAstro(x::Real, y::Real) = PhysicalVector2D(x*u"kpc", y*u"kpc")
PositionAstro(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"kpc", y*u"kpc", z*u"kpc")

############      Position       ###########

Velocity() = PhysicalVector3D(0.0u"m/s", 0.0u"m/s", 0.0u"m/s")
Velocity(u::Units=u"m/s") = PhysicalVector3D(0.0u, 0.0u, 0.0u)
Velocity(x::Real, y::Real, u::Units=u"m/s") = PhysicalVector2D(x*u, y*u)
Velocity(x::Real, y::Real, z::Real, u::Units=u"m/s") = PhysicalVector3D(x*u, y*u, z*u)
Velocity(x::Real, y::Real) = PhysicalVector2D(x*u"m/s", y*u"m/s")
Velocity(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"m/s", y*u"m/s", z*u"m/s")

VelocityAstro() = PhysicalVector3D(0.0u"kpc/Gyr", 0.0u"kpc/Gyr", 0.0u"kpc/Gyr")
VelocityAstro(u::Units=u"kpc/Gyr") = PhysicalVector3D(0.0u, 0.0u, 0.0u)
VelocityAstro(x::Real, y::Real, u::Units=u"kpc/Gyr") = PhysicalVector2D(x*u, y*u)
VelocityAstro(x::Real, y::Real, z::Real, u::Units=u"kpc/Gyr") = PhysicalVector3D(x*u, y*u, z*u)
VelocityAstro(x::Real, y::Real) = PhysicalVector2D(x*u"kpc/Gyr", y*u"kpc/Gyr")
VelocityAstro(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"kpc/Gyr", y*u"kpc/Gyr", z*u"kpc/Gyr")

############      Acceleration       ###########

Acceleration() = PhysicalVector3D(0.0u"m/s^2", 0.0u"m/s^2", 0.0u"m/s^2")
Acceleration(u::Units=u"m/s^2") = PhysicalVector3D(0.0u, 0.0u, 0.0u)
Acceleration(x::Real, y::Real, u::Units=u"m/s^2") = PhysicalVector2D(x*u, y*u)
Acceleration(x::Real, y::Real, z::Real, u::Units=u"m/s^2") = PhysicalVector3D(x*u, y*u, z*u)
Acceleration(x::Real, y::Real) = PhysicalVector2D(x*u"m/s^2", y*u"m/s^2")
Acceleration(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"m/s^2", y*u"m/s^2", z*u"m/s^2")

AccelerationAstro() = PhysicalVector3D(0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2")
AccelerationAstro(u::Units=u"kpc/Gyr^2") = PhysicalVector3D(0.0u, 0.0u, 0.0u)
AccelerationAstro(x::Real, y::Real, u::Units=u"kpc/Gyr^2") = PhysicalVector2D(x*u, y*u)
AccelerationAstro(x::Real, y::Real, z::Real, u::Units=u"kpc/Gyr^2") = PhysicalVector3D(x*u, y*u, z*u)
AccelerationAstro(x::Real, y::Real) = PhysicalVector2D(x*u"kpc/Gyr^2", y*u"kpc/Gyr^2")
AccelerationAstro(x::Real, y::Real, z::Real) = PhysicalVector3D(x*u"kpc/Gyr^2", y*u"kpc/Gyr^2", z*u"kpc/Gyr^2")

############      Basic mathematics       ###########
@inline length(p::T) where T <: AbstractPoint = 1
@inline iterate(p::T) where T <: AbstractPoint = (p,nothing)
@inline iterate(p::T,st) where T <: AbstractPoint = nothing
@inline real(p::T) where T <: AbstractPoint = p

# 2D
@inline +(p1::PhysicalVector2D, p2::PhysicalVector2D) = PhysicalVector2D(p1.x+p2.x, p1.y+p2.y)
@inline -(p1::PhysicalVector2D, p2::PhysicalVector2D) = PhysicalVector2D(p1.x-p2.x, p1.y-p2.y)
@inline *(p1::PhysicalVector2D, p2::PhysicalVector2D) = p1.x*p2.x + p1.y*p2.y
@inline +(p1::PhysicalVector2D, a::Quantity) = PhysicalVector2D(p1.x+a, p1.y+a)
@inline -(p1::PhysicalVector2D, a::Quantity) = PhysicalVector2D(p1.x-a, p1.y-a)
@inline *(p1::PhysicalVector2D, a::Quantity) = PhysicalVector2D(p1.x*a, p1.y*a)
@inline /(p1::PhysicalVector2D, a::Quantity) = PhysicalVector2D(p1.x/a, p1.y/a)
@inline *(p1::PhysicalVector2D, a::Real) = PhysicalVector2D(p1.x*a, p1.y*a)
@inline /(p1::PhysicalVector2D, a::Real) = PhysicalVector2D(p1.x/a, p1.y/a)
@inline +(a::Quantity, p::PhysicalVector2D) = PhysicalVector2D(a+p.x, a+p.y)
@inline -(a::Quantity, p::PhysicalVector2D) = PhysicalVector2D(a-p.x, a-p.y)
@inline *(a::Quantity, p::PhysicalVector2D) = PhysicalVector2D(a*p.x, a*p.y)
@inline /(a::Quantity, p::PhysicalVector2D) = PhysicalVector2D(a/p.x, a/p.y)
@inline *(a::Real, p::PhysicalVector2D) = PhysicalVector2D(a*p.x, a*p.y)
@inline norm(p::AbstractPoint2D) = sqrt(p.x^2 + p.y^2)
@inline dot(p1::PhysicalVector2D, p2::PhysicalVector2D) = p1.x*p2.x + p1.y*p2.y
@inline zero(p::PhysicalVector2D) = PhysicalVector2D(p.x*0.0, p.y*0.0)

# 3D
@inline +(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z)
@inline -(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
@inline *(p1::PhysicalVector3D, p2::PhysicalVector3D) = p1.x*p2.x + p1.y*p2.y + p1.z*p2.z
@inline *(p1::PhysicalVector3D, a::Real) = PhysicalVector3D(p1.x*a, p1.y*a, p1.z*a)
@inline +(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x+a, p1.y+a, p1.z+a)
@inline -(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x-a, p1.y-a, p1.z-a)
@inline *(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x*a, p1.y*a, p1.z*a)
@inline /(p1::PhysicalVector3D, a::Quantity) = PhysicalVector3D(p1.x/a, p1.y/a, p1.z/a)
@inline *(p1::PhysicalVector3D, a::Real) = PhysicalVector3D(p1.x*a, p1.y*a, p1.z*a)
@inline /(p1::PhysicalVector3D, a::Real) = PhysicalVector3D(p1.x/a, p1.y/a, p1.z/a)
@inline +(a::Quantity, p::PhysicalVector3D) = PhysicalVector3D(a+p.x, a+p.y, a+p.z)
@inline -(a::Quantity, p::PhysicalVector3D) = PhysicalVector3D(a-p.x, a-p.y, a-p.z)
@inline *(a::Quantity, p::PhysicalVector3D) = PhysicalVector3D(a*p.x, a*p.y, a*p.z)
@inline /(a::Quantity, p::PhysicalVector3D) = PhysicalVector3D(a/p.x, a/p.y, a/p.z)
@inline *(a::Real, p::PhysicalVector3D) = PhysicalVector3D(a*p.x, a*p.y, a*p.z)
@inline norm(p::AbstractPoint3D) = sqrt(p.x^2 + p.y^2 + p.z^2)
@inline dot(p1::PhysicalVector3D, p2::PhysicalVector3D) = p1.x*p2.x + p1.y*p2.y + p1.z*p2.z
@inline zero(p::PhysicalVector3D) = PhysicalVector3D(p.x*0.0, p.y*0.0, p.z*0.0)
@inline cross(p1::PhysicalVector3D, p2::PhysicalVector3D) = PhysicalVector3D(p1.y*p2.z-p1.z*p2.y, p1.z*p2.x-p1.x*p2.z, p1.x*p2.y-p1.y*p2.x)

############      Linear Algebra       ###########
@inline normalize(p::PhysicalVector2D) = (n = ustrip(norm(p)); return PhysicalVector(p.x/n, p.y/n))
@inline normalize(p::PhysicalVector3D) = (n = ustrip(norm(p)); return PhysicalVector(p.x/n, p.y/n, p.z/n))

@inline rotate_z(p::PhysicalVector2D, theta::Float64) = PhysicalVector2D(p.x*cos(theta)-p.y*sin(theta), p.x*sin(theta)+p.y*cos(theta))
@inline rotate(p::PhysicalVector2D, theta::Float64) = rotate_z(p, theta)

@inline rotate_x(p::PhysicalVector3D, theta::Float64) = PhysicalVector3D(p.x, p.y*cos(theta)-p.z*sin(theta), p.y*sin(theta)+p.z*cos(theta))
@inline rotate_y(p::PhysicalVector3D, theta::Float64) = PhysicalVector3D(p.x*cos(theta)+p.z*sin(theta), p.y, -p.x*sin(theta)+p.z*cos(theta))
@inline rotate_z(p::PhysicalVector3D, theta::Float64) = PhysicalVector3D(p.x*cos(theta)-p.y*sin(theta), p.x*sin(theta)+p.y*cos(theta), p.z)

"Computes the mean vector of an array of vectors"
function mean(a::Array{PhysicalVector3D})
    len = length(a)
    p = a[1]
    for i in 2:len
        @inbounds p += a[i]
    end
    return p/len
end

function mean(a::Array{PhysicalVector2D})
    len = length(a)
    p = a[1]
    for i in 2:len
        @inbounds p += a[i]
    end
    return p/len
end

"Converts Number array to PhysicalVector array"
function pconvert(a::Array{Float64,1}, u::Units)
    if length(a) == 3
        return PhysicalVector3D(a[1], a[2], a[3], u)
    elseif length(a) == 2
        return PhysicalVector2D(a[1], a[2], u)
    else
        error("Not supported dimension!")
    end
end

function pconvert(a::Array{Float64,1}, u::Units)
    if length(a) == 3
        return PhysicalVector3D(a[1], a[2], a[3], u)
    elseif length(a) == 2
        return PhysicalVector2D(a[1], a[2], u)
    else
        error("Not supported dimension!")
    end
end

function pconvert(a::Array{Float64,2}, u::Units)
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
            @inbounds append!(p, PhysicalVector2D(a[1,i], a[2,i], u))
        end
        return p
    else
        error("Not supported dimension!")
    end
end

############      Physical Particles       ###########
@enum ParticleType begin
    gas = 1
    halo = 2
    disk = 3
    bulge = 4
    star = 5
    blackhole = 6
end

abstract type AbstractParticle end
abstract type AbstractParticle2D <: AbstractParticle end
abstract type AbstractParticle3D <: AbstractParticle end
mutable struct PhysicalParticle2D <: AbstractParticle2D
    Pos::PhysicalVector2D
    Vel::PhysicalVector2D
    Acc::PhysicalVector2D
    Mass::Quantity
    ID::Int64
    Type::ParticleType
end
PhysicalParticle2D() = PhysicalParticle2D(PositionAstro(0.0,0.0), VelocityAstro(0.0,0.0), AccelerationAstro(0.0,0.0),
                                            0.0u"Msun", 0, ParticleType(5))

mutable struct PhysicalParticle3D <: AbstractParticle3D
    Pos::PhysicalVector3D
    Vel::PhysicalVector3D
    Acc::PhysicalVector3D
    Mass::Quantity
    ID::Int64
    Type::ParticleType
end
PhysicalParticle3D() = PhysicalParticle3D(PositionAstro(0.0,0.0,0.0), VelocityAstro(0.0,0.0,0.0), AccelerationAstro(0.0,0.0,0.0),
                                            0.0u"Msun", 0, ParticleType(5))

PhysicalParticle = PhysicalParticle3D

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
GasParticle2D() = GasParticle2D(PositionAstro(0.0,0.0), VelocityAstro(0.0,0.0), AccelerationAstro(0.0,0.0),
                                0.0u"Msun", 0,
                                0.0u"J/K", 0.0u"Msun/kpc^2", 0.0u"kpc",
                                0.0, 0.0, 0,
                                VelocityAstro(0.0,0.0), 0.0u"Gyr^-1", 0.0u"Gyr^-1", 0.0,
                                0.0u"N/m", 0.0u"J/K/s", 0.0u"kpc/s")

mutable struct GasData2D
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

GasData2D() = GasData2D(0.0u"J/K", 0.0u"Msun/kpc^2", 0.0u"kpc",
                        0.0, 0.0, 0,
                        VelocityAstro(0.0,0.0), 0.0u"Gyr^-1", 0.0u"Gyr^-1", 0.0,
                        0.0u"N/m", 0.0u"J/K/s", 0.0u"kpc/s")

mutable struct GasParticle3D
    Pos::PhysicalVector3D
    Vel::PhysicalVector3D
    Acc::PhysicalVector3D
    Mass::Quantity
    ID::Int64

    Entropy::Quantity
    Density::Quantity
    Hsml::Quantity

    Left::Float64
    Right::Float64
    NumNgbFound::Int64

    RotVel::PhysicalVector3D
    DivVel::Quantity
    CurlVel::Quantity
    dHsmlRho::Float64

    Pressure::Quantity
    DtEntropy::Quantity
    MaxSignalVel::Quantity
end
GasParticle3D() = GasParticle3D(PositionAstro(0.0,0.0,0.0), VelocityAstro(0.0,0.0,0.0), AccelerationAstro(0.0,0.0,0.0),
                                0.0u"Msun", 0,
                                0.0u"J/K", 0.0u"Msun/kpc^3", 0.0u"kpc",
                                0.0, 0.0, 0,
                                VelocityAstro(0.0,0.0,0.0), 0.0u"Gyr^-1", 0.0u"Gyr^-1", 0.0,
                                0.0u"N/m^2", 0.0u"J/K/s", 0.0u"kpc/s")

GasParticle = GasParticle3D

mutable struct GasData3D
    Entropy::Quantity
    Density::Quantity
    Hsml::Quantity

    Left::Float64
    Right::Float64
    NumNgbFound::Int64

    RotVel::PhysicalVector3D
    DivVel::Quantity
    CurlVel::Quantity
    dHsmlRho::Float64

    Pressure::Quantity
    DtEntropy::Quantity
    MaxSignalVel::Quantity
end
GasData3D() = GasData3D(0.0u"J/K", 0.0u"Msun/kpc^3", 0.0u"kpc",
                        0.0, 0.0, 0,
                        VelocityAstro(0.0,0.0,0.0), 0.0u"Gyr^-1", 0.0u"Gyr^-1", 0.0,
                        0.0u"N/m^2", 0.0u"J/K/s", 0.0u"kpc/s")

GasData = GasData3D

############      Constants      ###########
@constant(H, "Hubble constant", 74.03, BigFloat(74.03),
            u"km/s/Mpc", 1.42, BigFloat(1.42), "Hubble Space Telescope 2019-03-18")
@constant(ACC0, "Modified gravitational acceleration", 1.2e-8, BigFloat(1.2e-8),
            u"cm/s^2", 0.0, BigFloat(0.0), "Milgrom 1983")
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
PhysicalConstant() = PhysicalConstant(CODATA2014.c,
                                      CODATA2014.G,
                                      CODATA2014.h,
                                      CODATA2014.e,
                                      CODATA2014.m_e,
                                      CODATA2014.m_n,
                                      CODATA2014.m_p,
                                      CODATA2014.σ,
                                      H,
                                      ACC0)
Constants = PhysicalConstant()

mutable struct Header_Gadget2 # Refer to Gadget2 manual for more information
    npart::Array{Int32,1} # gas, halo, disk, Bulge, star, blackholw
    mass::Array{Float64,1}

    time::Float64
    redshift::Float64

    flag_sfr::Int32
    flag_feedback::Int32

    npartTotal::Array{UInt32,1}

    flag_cooling::Int32

    num_files::Int32

    BoxSize::Float64
    Omega0::Float64
    OmegaLambda::Float64
    HubbleParam::Float64

    flag_stellarage::Int32
    flag_metals::Int32

    npartTotalHighWord::Array{UInt32,1}

    flag_entropy_instead_u::Int32
    # fill 60 char

    # Some
end # Header_Gadget2

Header_Gadget2() = Header_Gadget2([0,0,0,0,0,0],
                                    [0.0,0.0,0.0,0.0,0.0,0.0],
                                    0.0, 0.0, 0, 0,
                                    [0,0,0,0,0,0],
                                    0, 1, 0.0, 0.3, 0.7, 0.71, 0, 0,
                                    [0,0,0,0,0,0], 0)

mutable struct TreeNode
    ID::Int64
    Father::Int64
    DaughterID::Array{Int64,1}
    Center::Array{Float64,1}
    SideLength::Float64
    MaxSoft::Float64
    SparseDaughterID::Array{Int64,1} # Walk in sparse tree to improve performance
    IsLeaf::Bool
    ParticleID::Int64 # Refers to the particle on this leaf.
                      # One leaf can only take one particle. Set 0 if none or more than 1
end # TreeNode

mutable struct PhysicalTreeNode
    ID::Int64
    Father::Int64
    DaughterID::Array{Int64,1}
    Center::AbstractPoint
    SideLength::Quantity
    Mass::Quantity
    MassCenter::AbstractPoint
    MaxSoft::Quantity
    SparseDaughterID::Array{Int64,1} # Walk in sparse tree to improve performance
    PDM_Mass::Quantity
    PDM_MassCenter::AbstractPoint
    IsLeaf::Bool
    ParticleID::Int64 # Refers to the particle on this leaf.
                      # One leaf can only take one particle. Set 0 if none or more than 1
end # PhysicalTreeNode

include("io/write_ascii.jl")
include("io/write_gadget.jl")

include("io/read_ascii.jl")
include("io/read_gadget.jl")

include("numerics/basic.jl")
include("numerics/center.jl")
include("numerics/kdtree.jl")
include("numerics/peano.jl")

end  # module PhysicalParticles
