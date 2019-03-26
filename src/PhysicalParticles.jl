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

"Returns the minimum x value of an array of points"
function min_x(a::Array{T,1}) where T <: AbstractPoint
    min = a[1].x
    for p in a
        if min > p.x
            min = p.x
        end
    end
    return min
end

"Returns the minimum y value of an array of points"
function min_y(a::Array{T,1}) where T <: AbstractPoint
    min = a[1].y
    for p in a
        if min > p.y
            min = p.y
        end
    end
    return min
end

"Returns the minimum z value of an array of points"
function min_z(a::Array{T,1}) where T <: AbstractPoint3D
    min = a[1].z
    for p in a
        if min > p.z
            min = p.z
        end
    end
    return min
end

"Returns the maximum x value of an array of points"
function max_x(a::Array{T,1}) where T <: AbstractPoint
    max = a[1].x
    for p in a
        if max < p.x
            max = p.x
        end
    end
    return max
end

"Returns the maximum y value of an array of points"
function max_y(a::Array{T,1}) where T <: AbstractPoint
    max = a[1].y
    for p in a
        if max < p.y
            max = p.y
        end
    end
    return max
end

"Returns the maximum z value of an array of points"
function max_z(a::Array{T,1}) where T <: AbstractPoint3D
    max = a[1].z
    for p in a
        if max < p.z
            max = p.z
        end
    end
    return max
end

"Returns x center of the box"
function center_x(a::Array{T,1}) where T <: AbstractPoint
    left = min_x(a)
    right = max_x(a)
    return (left + right) / 2.0
end

"Returns y center of the box"
function center_y(a::Array{T,1}) where T <: AbstractPoint
    left = min_y(a)
    right = max_y(a)
    return (left + right) / 2.0
end

"Returns z center of the box"
function center_z(a::Array{T,1}) where T <: AbstractPoint3D
    left = min_z(a)
    right = max_z(a)
    return (left + right) / 2.0
end

"Returns center of the box"
function center(a::Array{T,1}) where T <: AbstractPoint2D
    x = center_x(a)
    y = center_y(a)
    return typeof(a[1])(x,y)
end

function center(a::Array{T,1}) where T <: AbstractPoint3D
    x = center_x(a)
    y = center_y(a)
    z = center_z(a)
    return typeof(a[1])(x,y,z)
end

#### Extent
# Seperate normal points and physical vectors for type stabilization
struct Extent2D
    xMin::Float64
    xMax::Float64
    yMin::Float64
    yMax::Float64
    SideLength::Float64
    Center::Point2D
end

struct Extent3D
    xMin::Float64
    xMax::Float64
    yMin::Float64
    yMax::Float64
    zMin::Float64
    zMax::Float64
    SideLength::Float64
    Center::Point3D
end

"Returns extent of normal particles"
Extent(a::Array{Point2D}) = (xMin=min_x(a); xMax=max_x(a); yMin=min_y(a); yMax=max_y(a);
                              len=max(xMax-xMin, yMax-yMin);
                              Center=Point2D(0.5(xMax-xMin), 0.5(yMax-yMin));
                              return Extent2D(xMin,xMax,yMin,yMax,len,Center))

Extent(a::Array{Point3D}) = (xMin=min_x(a); xMax=max_x(a); yMin=min_y(a); yMax=max_y(a); zMin=min_z(a); zMax=max_z(a);
                              len=max(xMax-xMin, yMax-yMin, zMax-zMin);
                              Center=Point3D(0.5(xMax-xMin), 0.5(yMax-yMin), 0.5(zMax-zMin));
                              return Extent3D(xMin,xMax,yMin,yMax,zMin,zMax,len,Center))



struct PhysicalExtent2D
    xMin::Quantity
    xMax::Quantity
    yMin::Quantity
    yMax::Quantity
    SideLength::Quantity
    Center::PhysicalVector2D
end

struct PhysicalExtent3D
    xMin::Quantity
    xMax::Quantity
    yMin::Quantity
    yMax::Quantity
    zMin::Quantity
    zMax::Quantity
    SideLength::Quantity
    Center::PhysicalVector3D
end

"Returns extent of normal particles"
PhysicalExtent(a::Array{Point2D}) = (xMin=min_x(a); xMax=max_x(a); yMin=min_y(a); yMax=max_y(a);
                                        len=max(xMax-xMin, yMax-yMin);
                                        Center=PhysicalVector2D(0.5(xMax-xMin), 0.5(yMax-yMin));
                                        return PhysicalExtent2D(xMin,xMax,yMin,yMax,len,Center))

PhysicalExtent(a::Array{Point3D}) = (xMin=min_x(a); xMax=max_x(a); yMin=min_y(a); yMax=max_y(a); zMin=min_z(a); zMax=max_z(a);
                                        len=max(xMax-xMin, yMax-yMin, zMax-zMin);
                                        Center=PhysicalVector3D(0.5(xMax-xMin), 0.5(yMax-yMin), 0.5(zMax-zMin));
                                        return PhysicalExtent3D(xMin,xMax,yMin,yMax,zMin,zMax,len,Center))

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

mutable struct PhysicalParticle3D <: AbstractParticle3D
    Pos::PhysicalVector3D
    Vel::PhysicalVector3D
    Acc::PhysicalVector3D
    Mass::Quantity
    ID::Int64
    Type::Int64
end

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

    RotVel::PhysicalVector2D
    DivVel::Quantity
    CurlVel::Quantity
    dHsmlRho::Quantity

    Pressure::Quantity
    DtEntropy::Quantity
    MaxSignalVel::Quantity
end

mutable struct GasParticle3D <: AbstractParticle3D
    Pos::PhysicalVector3D
    Vel::PhysicalVector3D
    Acc::PhysicalVector3D
    Mass::Quantity
    ID::Int64

    Entropy::Quantity
    Density::Quantity
    Hsml::Quantity

    RotVel::PhysicalVector3D
    DivVel::Quantity
    CurlVel::Quantity
    dHsmlRho::Quantity

    Pressure::Quantity
    DtEntropy::Quantity
    MaxSignalVel::Quantity
end

############      Constants      ###########

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
    PhysicalConstant() = PhysicalConstant(CODATA2014.c,
                                          CODATA2014.G,
                                          CODATA2014.h,
                                          CODATA2014.e,
                                          CODATA2014.m_e,
                                          CODATA2014.m_n,
                                          CODATA2014.m_p,
                                          CODATA2014.σ,
                                          Constant(H, "Hubble constant", 74.03, BigFloat(74.03),
                                                    u"km/s/Mpc", 1.42, BigFloat(1.42), "Hubble Space Telescope 2019-03-18") # Wiki, 2019-03-18
                                          )
end

############      Peano-Hilbert       ###########
# Copied from GeometicalPredicates.jl and referred to Gadget2

# number of bits to use per dimension in calculating the peano-key
const peano_2D_bits = 31
const peano_3D_bits = 21

# implementing 2D scale dependednt Peano-Hilbert indexing

 _extract_peano_bin_num(nbins::Int64, n::Float64) = trunc(Integer, (n-1)*nbins )

# calculate peano key for given point
function peanokey(p::PhysicalVector2D, bits::Int64=peano_2D_bits)
    n = 1 << bits
    s = n >> 1; d = 0
    x = _extract_peano_bin_num(n, ustrip(Float64, u, p.x))
    y = _extract_peano_bin_num(n, ustrip(Float64, u, p.y))
    while true
        rx = (x & s) > 0
        ry = (y & s) > 0
        d += s * s * xor(3 * rx, ry)
        s = s >> 1
        (s == 0) && break
        if ry == 0
            if rx == 1
                x = n - 1 - x;
                y = n - 1 - y;
            end
            x, y = y, x
        end
    end
end

# implementing 3D scaleful Peano-Hilbert indexing
const quadrants_arr = [
  0, 7, 1, 6, 3, 4, 2, 5,
  7, 4, 6, 5, 0, 3, 1, 2,
  4, 3, 5, 2, 7, 0, 6, 1,
  3, 0, 2, 1, 4, 7, 5, 6,
  1, 0, 6, 7, 2, 3, 5, 4,
  0, 3, 7, 4, 1, 2, 6, 5,
  3, 2, 4, 5, 0, 1, 7, 6,
  2, 1, 5, 6, 3, 0, 4, 7,
  6, 1, 7, 0, 5, 2, 4, 3,
  1, 2, 0, 3, 6, 5, 7, 4,
  2, 5, 3, 4, 1, 6, 0, 7,
  5, 6, 4, 7, 2, 1, 3, 0,
  7, 6, 0, 1, 4, 5, 3, 2,
  6, 5, 1, 2, 7, 4, 0, 3,
  5, 4, 2, 3, 6, 7, 1, 0,
  4, 7, 3, 0, 5, 6, 2, 1,
  6, 7, 5, 4, 1, 0, 2, 3,
  7, 0, 4, 3, 6, 1, 5, 2,
  0, 1, 3, 2, 7, 6, 4, 5,
  1, 6, 2, 5, 0, 7, 3, 4,
  2, 3, 1, 0, 5, 4, 6, 7,
  3, 4, 0, 7, 2, 5, 1, 6,
  4, 5, 7, 6, 3, 2, 0, 1,
  5, 2, 6, 1, 4, 3, 7, 0]
quadrants(a::Int64, b::Int64, c::Int64, d::Int64) = (@inbounds x = quadrants_arr[1+a<<3+b<<2+c<<1+d]; x)
rotxmap_table = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22]
rotymap_table = [1, 2, 3, 0, 16, 17, 18, 19, 11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7]
rotx_table = [3, 0, 0, 2, 2, 0, 0, 1]
roty_table = [0, 1, 1, 2, 2, 3, 3, 0]
sense_table = [-1, -1, -1, +1, +1, -1, -1, -1]

function peanokey(p::PhysicalVector3D, u::Units, bits::Int64=peano_3D_bits)
    n = 1 << bits
    x = _extract_peano_bin_num(n, ustrip(Float64, u, p.x))
    y = _extract_peano_bin_num(n, ustrip(Float64, u, p.y))
    z = _extract_peano_bin_num(n, ustrip(Float64, u, p.z))
    mask = 1 << (bits - 1)
    key = 0
    rotation = 0
    sense = 1
    for i in 1:bits
        bitx = (x & mask > 0) ? 1 : 0
        bity = (y & mask > 0) ? 1 : 0
        bitz = (z & mask > 0) ? 1 : 0

        quad = quadrants(rotation, bitx, bity, bitz)

        key <<= 3
        key += sense == 1 ? quad : 7-quad

        @inbounds rotx = rotx_table[quad+1]
        @inbounds roty = roty_table[quad+1]
        @inbounds sense *= sense_table[quad+1]

        while rotx > 0
            @inbounds rotation = rotxmap_table[rotation+1]
            rotx -= 1
        end

        while(roty > 0)
            @inbounds rotation = rotymap_table[rotation+1]
            roty -= 1
        end
        mask >>= 1
    end

    key
end

# implementing scale-free Hilbert ordering. Real all about it here:
# http://doc.cgal.org/latest/Spatial_sorting/index.html

abstract type AbstractCoordinate end
mutable struct CoordinateX <: AbstractCoordinate end
mutable struct CoordinateY <: AbstractCoordinate end
mutable struct CoordinateZ <: AbstractCoordinate end
const coordinatex = CoordinateX()
const coordinatey = CoordinateY()
const coordinatez = CoordinateZ()
next2d(::CoordinateX) = coordinatey
next2d(::CoordinateY) = coordinatex
next3d(::CoordinateX) = coordinatey
next3d(::CoordinateY) = coordinatez
next3d(::CoordinateZ) = coordinatex
nextnext3d(::CoordinateX) = coordinatez
nextnext3d(::CoordinateY) = coordinatex
nextnext3d(::CoordinateZ) = coordinatey

abstract type AbstractDirection end
mutable struct Forward <: AbstractDirection end
mutable struct Backward <: AbstractDirection end
const forward = Forward()
const backward = Backward()
Base.:!(::Forward) = backward
Base.:!(::Backward) = forward

compare(::Forward, ::CoordinateX, p1::AbstractPoint, p2::AbstractPoint) = getx(p1) < getx(p2)
compare(::Backward, ::CoordinateX, p1::AbstractPoint, p2::AbstractPoint) = getx(p1) > getx(p2)
compare(::Forward, ::CoordinateY, p1::AbstractPoint, p2::AbstractPoint) = gety(p1) < gety(p2)
compare(::Backward, ::CoordinateY, p1::AbstractPoint, p2::AbstractPoint) = gety(p1) > gety(p2)
compare(::Forward, ::CoordinateZ, p1::AbstractPoint, p2::AbstractPoint) = getz(p1) < getz(p2)
compare(::Backward, ::CoordinateZ, p1::AbstractPoint, p2::AbstractPoint) = getz(p1) > getz(p2)

function select!(direction::AbstractDirection, coordinate::AbstractCoordinate, v::Array{T,1}, k::Int, lo::Int, hi::Int) where T<:AbstractPoint
    lo <= k <= hi || error("select index $k is out of range $lo:$hi")
    @inbounds while lo < hi
        if hi-lo == 1
            if compare(direction, coordinate, v[hi], v[lo])
                v[lo], v[hi] = v[hi], v[lo]
            end
            return v[k]
        end
        pivot = v[(lo+hi)>>>1]
        i, j = lo, hi
        while true
            while compare(direction, coordinate, v[i], pivot); i += 1; end
            while compare(direction, coordinate, pivot, v[j]); j -= 1; end
            i <= j || break
            v[i], v[j] = v[j], v[i]
            i += 1; j -= 1
        end
        if k <= j
            hi = j
        elseif i <= k
            lo = i
        else
            return pivot
        end
    end
    return v[lo]
end

function hilbertsort!(directionx::AbstractDirection, directiony::AbstractDirection, coordinate::AbstractCoordinate, a::Array{T,1}, lo::Int64, hi::Int64, lim::Int64=4) where T<:AbstractPoint2D
    hi-lo <= lim && return a

    i2 = (lo+hi)>>>1
    i1 = (lo+i2)>>>1
    i3 = (i2+hi)>>>1

    select!(directionx, coordinate, a, i2, lo, hi)
    select!(directiony, next2d(coordinate), a, i1, lo, i2)
    select!(!directiony, next2d(coordinate), a, i3, i2, hi)

    hilbertsort!(directiony, directionx, next2d(coordinate), a, lo, i1, lim)
    hilbertsort!(directionx, directiony, coordinate, a, i1, i2, lim)
    hilbertsort!(directionx, directiony, coordinate, a, i2, i3, lim)
    hilbertsort!(!directiony, !directionx, next2d(coordinate), a, i3, hi, lim)

    return a
end

function hilbertsort!(directionx::AbstractDirection, directiony::AbstractDirection, directionz::AbstractDirection, coordinate::AbstractCoordinate, a::Array{T,1}, lo::Int64, hi::Int64, lim::Int64=8) where T<:AbstractPoint3D
    hi-lo <= lim && return a

    i4 = (lo+hi)>>>1
    i2 = (lo+i4)>>>1
    i1 = (lo+i2)>>>1
    i3 = (i2+i4)>>>1
    i6 = (i4+hi)>>>1
    i5 = (i4+i6)>>>1
    i7 = (i6+hi)>>>1

    select!(directionx, coordinate, a, i4, lo, hi)
    select!(directiony, next3d(coordinate), a, i2, lo, i4)
    select!(directionz, nextnext3d(coordinate), a, i1, lo, i2)
    select!(!directionz, nextnext3d(coordinate), a, i3, i2, i4)
    select!(!directiony, next3d(coordinate), a, i6, i4, hi)
    select!(directionz, nextnext3d(coordinate), a, i5, i4, i6)
    select!(!directionz, nextnext3d(coordinate), a, i7, i6, hi)

    hilbertsort!( directionz,  directionx,  directiony, nextnext3d(coordinate), a, lo, i1, lim)
    hilbertsort!( directiony,  directionz,  directionx, next3d(coordinate),     a, i1, i2, lim)
    hilbertsort!( directiony,  directionz,  directionx, next3d(coordinate),     a, i2, i3, lim)
    hilbertsort!( directionx, !directiony, !directionz, coordinate,             a, i3, i4, lim)
    hilbertsort!( directionx, !directiony, !directionz, coordinate,             a, i4, i5, lim)
    hilbertsort!(!directiony,  directionz, !directionx, next3d(coordinate),     a, i5, i6, lim)
    hilbertsort!(!directiony,  directionz, !directionx, next3d(coordinate),     a, i6, i7, lim)
    hilbertsort!(!directionz, !directionx,  directiony, nextnext3d(coordinate), a, i7, hi, lim)

    return a
end

hilbertsort!(a::Array{T,1}) where {T<:AbstractPoint2D} = hilbertsort!(backward, backward, coordinatey, a, 1, length(a))
hilbertsort!(a::Array{T,1}, lo::Int64, hi::Int64, lim::Int64) where {T<:AbstractPoint2D} = hilbertsort!(backward, backward, coordinatey, a, lo, hi, lim)
hilbertsort!(a::Array{T,1}) where {T<:AbstractPoint3D} = hilbertsort!(backward, backward, backward, coordinatez, a, 1, length(a))
hilbertsort!(a::Array{T,1}, lo::Int64, hi::Int64, lim::Int64) where {T<:AbstractPoint3D} = hilbertsort!(backward, backward, backward, coordinatey, a, lo, hi, lim)

# multi-scale sort. Read all about it here:
# http://doc.cgal.org/latest/Spatial_sorting/classCGAL_1_1Multiscale__sort.html
function _mssort!(a::Array{T,1}, lim_ms::Int64, lim_hl::Int64, rat::Float64) where T<:AbstractPoint
    hi = length(a)
    lo = 1
    while true
        lo = hi - round(Int, (1-rat)*hi)
        hi-lo <= lim_ms && return a
        hilbertsort!(a, lo, hi, lim_hl)
        hi = lo-1
    end
    return a
end

# Utility methods, setting some different defaults for 2D and 3D. These are exported
mssort!(a::Array{T,1}; lim_ms::Int64=16, lim_hl::Int64=4, rat::Float64=0.25) where {T<:AbstractPoint2D} =
    _mssort!(a, lim_ms, lim_hl, rat)
mssort!(a::Array{T,1}; lim_ms::Int64=64, lim_hl::Int64=8, rat::Float64=0.125) where {T<:AbstractPoint3D} =
    _mssort!(a, lim_ms, lim_hl, rat)


### Normal Points

function peanokey(p::Point2D, bits::Int64=peano_2D_bits)
    n = 1 << bits
    s = n >> 1; d = 0
    x = _extract_peano_bin_num(n, getx(p))
    y = _extract_peano_bin_num(n, gety(p))
    while true
        rx = (x & s) > 0
        ry = (y & s) > 0
        d += s * s * xor(3 * rx, ry)
        s = s >> 1
        (s == 0) && break
        if ry == 0
            if rx == 1
                x = n - 1 - x;
                y = n - 1 - y;
            end
            x, y = y, x
        end
    end
    d
end

function peanokey(p::Point3D, bits::Int64=peano_3D_bits)
    n = 1 << bits
    x = _extract_peano_bin_num(n, getx(p))
    y = _extract_peano_bin_num(n, gety(p))
    z = _extract_peano_bin_num(n, getz(p))
    mask = 1 << (bits - 1)
    key = 0
    rotation = 0
    sense = 1
    for i in 1:bits
        bitx = (x & mask > 0) ? 1 : 0
        bity = (y & mask > 0) ? 1 : 0
        bitz = (z & mask > 0) ? 1 : 0

        quad = quadrants(rotation, bitx, bity, bitz)

        key <<= 3
        key += sense == 1 ? quad : 7-quad

        @inbounds rotx = rotx_table[quad+1]
        @inbounds roty = roty_table[quad+1]
        @inbounds sense *= sense_table[quad+1]

        while rotx > 0
            @inbounds rotation = rotxmap_table[rotation+1]
            rotx -= 1
        end

        while(roty > 0)
            @inbounds rotation = rotymap_table[rotation+1]
            roty -= 1
        end
        mask >>= 1
    end

    key
end
end  # module PhysicalParticles
