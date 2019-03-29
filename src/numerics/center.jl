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

PhysicalExtent(a::Array{PhysicalVector2D}) = (xMin=min_x(a); xMax=max_x(a); yMin=min_y(a); yMax=max_y(a);
                                              len=max(xMax-xMin, yMax-yMin);
                                              Center=PhysicalVector2D(0.5(xMax-xMin), 0.5(yMax-yMin));
                                              return PhysicalExtent2D(xMin,xMax,yMin,yMax,len,Center))

PhysicalExtent(a::Array{PhysicalVector3D}) = (xMin=min_x(a); xMax=max_x(a); yMin=min_y(a); yMax=max_y(a); zMin=min_z(a); zMax=max_z(a);
                                              len=max(xMax-xMin, yMax-yMin, zMax-zMin);
                                              Center=PhysicalVector3D(0.5(xMax-xMin), 0.5(yMax-yMin), 0.5(zMax-zMin));
                                              return PhysicalExtent3D(xMin,xMax,yMin,yMax,zMin,zMax,len,Center))


#####  Mass center  #####

function mass_center(a::Array{T}) where T <: AbstractParticle
    if length(a) == 1
        return a.Pos
    end
    sum_mass = a[1].Mass
    sum_center = a[1].Pos * a[1].Mass
    for i in 2:length(a)
        sum_center += a[i].Pos * a[i].Mass
    end

    return sum_center/sum_mass
end
