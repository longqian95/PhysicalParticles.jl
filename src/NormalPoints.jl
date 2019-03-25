# We support normal support for mathematical applications

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

#######   Basic mathematics   #######

# 2D
@inline +(p1::Point2D, p2::Point2D) = Point2D(p1.x+p2.x, p1.y+p2.y)
@inline -(p1::Point2D, p2::Point2D) = Point2D(p1.x-p2.x, p1.y-p2.y)
@inline *(p1::Point2D, p2::Point2D) = Point2D(p1.x*p2.x + p1.y*p2.y)
@inline +(p1::Point2D, a::Real) = Point2D(p1.x+a, p1.y+a)
@inline -(p1::Point2D, a::Real) = Point2D(p1.x-a, p1.y-a)
@inline *(p1::Point2D, a::Real) = Point2D(p1.x*a, p1.y*a)
@inline /(p1::Point2D, a::Real) = Point2D(p1.x/a, p1.y/a)
@inline +(a::Real, p1::Point2D) = Point2D(p1.x+a, p1.y+a)
@inline -(a::Real, p1::Point2D) = Point2D(p1.x-a, p1.y-a)
@inline *(a::Real, p1::Point2D) = Point2D(p1.x*a, p1.y*a)
@inline /(a::Real, p1::Point2D) = Point2D(p1.x/a, p1.y/a)
@inline dot(p1::Point2D, p2::Point2D) = Point2D(p1.x*p2.x + p1.y*p2.y)
@inline zero(p::Point2D) = Point2D(p.x*0.0, p.y*0.0)

# 3D
@inline +(p1::Point3D, p2::Point3D) = Point3D(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z)
@inline -(p1::Point3D, p2::Point3D) = Point3D(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
@inline *(p1::Point3D, p2::Point3D) = Point3D(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z)
@inline *(p1::Point3D, a::Real) = Point3D(p1.x*a, p1.y*a, p1.z*a)
@inline +(p1::Point3D, a::Real) = Point3D(p1.x+a, p1.y+a, p1.z+a)
@inline -(p1::Point3D, a::Real) = Point3D(p1.x-a, p1.y-a, p1.z-a)
@inline *(p1::Point3D, a::Real) = Point3D(p1.x*a, p1.y*a, p1.z*a)
@inline /(p1::Point3D, a::Real) = Point3D(p1.x/a, p1.y/a, p1.z/a)
@inline +(a::Real, p1::Point3D) = Point3D(p1.x+a, p1.y+a, p1.z+a)
@inline -(a::Real, p1::Point3D) = Point3D(p1.x-a, p1.y-a, p1.z-a)
@inline *(a::Real, p1::Point3D) = Point3D(p1.x*a, p1.y*a, p1.z*a)
@inline /(a::Real, p1::Point3D) = Point3D(p1.x/a, p1.y/a, p1.z/a)
@inline dot(p1::Point3D, p2::Point3D) = Point3D(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z)
@inline zero(p::Point3D) = Point3D(p.x*0.0, p.y*0.0, p.z*0.0)
@inline cross(p1::Point3D, p2::Point3D) = Point3D(p1.y*p2.z-p1.z*p2.y, p1.z*p2.x-p1.x*p2.z, p1.x*p2.y-p1.y*p2.x)

############      Linear Algebra       ###########

@inline norm(p::Point2D) = sqrt(p.x*p.x + p.y*p.y + p.z*p.z)
@inline normalize(p::Point2D) = (n = norm(p); return Point(p.x/n, p.y/n, p.z/n))
@inline normalize!(p::Point2D) = (n = norm(p); p.x /= n; p.y /= n; p.z /= n)

@inline rotate_z(p::Point2D, theta::Float64) = Point2D(p.x*cos(theta)-p.y*sin(theta), p.x*sin(theta)+p.y*cos(theta))
@inline rotate(p::Point2D, theta::Float64) = rotate_z(p, theta)

@inline rotate_x(p::Point3D, theta::Float64) = Point3D(p.x, p.y*cos(theta)-p.z*sin(theta), p.y*sin(theta)+p.z*cos(theta))
@inline rotate_y(p::Point3D, theta::Float64) = Point3D(p.x*cos(theta)+p.z*sin(theta), p.y, -p.x*sin(theta)+p.z*cos(theta))
@inline rotate_z(p::Point3D, theta::Float64) = Point3D(p.x*cos(theta)-p.y*sin(theta), p.x*sin(theta)+p.y*cos(theta))

"Computes the mean vector of an array of vectors"
function mean(a::Array{Point3D})
    len = length(a)
    p = Point3D()
    for i in 1:len
        p += a[i]
    end
    return p/len
end

function mean(a::Array{Point2D})
    len = length(a)
    p = Point3D()
    for i in 1:len
        p += a[i]
    end
    return p/len
end

"Converts Number array to Point array"
function pconvert(a::Array{Float64,1})
    if length(a) == 3
        return Point3D(a[1], a[2], a[3], u"m")
    elseif length(a) == 2
        return Point(a[1], a[2], u"m")
    else
        error("Not supported dimension!")
    end
end

function pconvert(a::Array{Float64,1}, u::Units)
    if length(a) == 3
        return Point3D(a[1], a[2], a[3], u)
    elseif length(a) == 2
        return Point(a[1], a[2], u)
    else
        error("Not supported dimension!")
    end
end

function pconvert(a::Array{Float64,2}, u::Units)
    row, col = size(a)
    if row == 3
        p = rand(Point3D,0)
        for i in 1:col
            @inbounds append!(p, Point3D(a[1,i], a[2,i], a[3,i], u))
        end
        return p
    elseif row == 2
        p = rand(Point2D,0)
        for i in 1:col
            @inbounds append!(p, Point3D(a[1,i], a[2,i], u))
        end
        return p
    else
        error("Not supported dimension!")
    end
end
