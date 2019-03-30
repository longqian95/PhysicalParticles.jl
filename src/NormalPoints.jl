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

@inline +(p1::Point2D, a::Real) = Point2D(p1.x+a, p1.y+a)
@inline -(p1::Point2D, a::Real) = Point2D(p1.x-a, p1.y-a)
@inline +(a::Real, p::Point2D) = Point2D(a+p.x, a+p.y)
@inline -(a::Real, p::Point2D) = Point2D(a-p.x, a-p.y)

@inline +(p1::Point3D, a::Real) = Point3D(p1.x+a, p1.y+a, p1.z+a)
@inline -(p1::Point3D, a::Real) = Point3D(p1.x-a, p1.y-a, p1.z-a)
@inline +(a::Real, p::Point3D) = Point3D(a+p.x, a+p.y, a+p.z)
@inline -(a::Real, p::Point3D) = Point3D(a-p.x, a-p.y, a-p.z)

############      Linear Algebra       ###########

@inline norm(p::Point2D) = sqrt(p.x*p.x + p.y*p.y)
@inline normalize(p::Point2D) = (n = norm(p); return Point(p.x/n, p.y/n))
@inline normalize(p::Point3D) = (n = norm(p); return Point(p.x/n, p.y/n, p.z/n))

@inline rotate_z(p::Point2D, theta::Float64) = Point2D(p.x*cos(theta)-p.y*sin(theta), p.x*sin(theta)+p.y*cos(theta))
@inline rotate(p::Point2D, theta::Float64) = rotate_z(p, theta)

@inline rotate_x(p::Point3D, theta::Float64) = Point3D(p.x, p.y*cos(theta)-p.z*sin(theta), p.y*sin(theta)+p.z*cos(theta))
@inline rotate_y(p::Point3D, theta::Float64) = Point3D(p.x*cos(theta)+p.z*sin(theta), p.y, -p.x*sin(theta)+p.z*cos(theta))
@inline rotate_z(p::Point3D, theta::Float64) = Point3D(p.x*cos(theta)-p.y*sin(theta), p.x*sin(theta)+p.y*cos(theta), p.z)

"Converts Number array to Point array"
function npconvert(a::Array{Float64,1})
    if length(a) == 3
        return Point3D(a[1], a[2], a[3])
    elseif length(a) == 2
        return Point2D(a[1], a[2])
    else
        error("Not supported dimension!")
    end
end

function npconvert(a::Array{Float64,2})
    row, col = size(a)
    if row == 3
        p = rand(Point3D,0)
        for i in 1:col
            @inbounds append!(p, Point3D(a[1,i], a[2,i], a[3,i]))
        end
        return p
    elseif row == 2
        p = rand(Point2D,0)
        for i in 1:col
            @inbounds append!(p, Point2D(a[1,i], a[2,i]))
        end
        return p
    else
        error("Not supported dimension!")
    end
end
