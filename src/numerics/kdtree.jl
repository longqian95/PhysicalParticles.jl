using NearestNeighbors

#####  Kd-tree setup  #####

function kdtree_setup(a::Array{Point2D,1})
    NumTotal = length(a)
    Pos = zeros(2,NumTotal)
    for i in 1:NumTotal
        @inbounds Pos[1,i] = a[i].x
        @inbounds Pos[2,i] = a[i].y
    end
    return kdtree = KDTree(Pos)
end

function kdtree_setup(a::Array{Point3D,1})
    NumTotal = length(a)
    Pos = zeros(3,NumTotal)
    for i in 1:NumTotal
        @inbounds Pos[1,i] = a[i].x
        @inbounds Pos[2,i] = a[i].y
        @inbounds Pos[3,i] = a[i].z
    end
    return kdtree = KDTree(Pos)
end

function kdtree_setup(a::Array{PhysicalVector2D,1}; mode = "Astro")
    NumTotal = length(a)
    Pos = zeros(2,NumTotal)
    if mode == "Astro"
        for i in 1:NumTotal
            @inbounds Pos[1,i] = ustrip(u"kpc", a[i].x)
            @inbounds Pos[2,i] = ustrip(u"kpc", a[i].y)
        end
    elseif mode == "SI"
        for i in 1:NumTotal
            @inbounds Pos[1,i] = ustrip(u"m", a[i].x)
            @inbounds Pos[2,i] = ustrip(u"m", a[i].y)
        end
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
    return kdtree = KDTree(Pos)
end

function kdtree_setup(a::Array{PhysicalVector3D,1}; mode = "Astro")
    NumTotal = length(a)
    Pos = zeros(3,NumTotal)
    if mode == "Astro"
        for i in 1:NumTotal
            @inbounds Pos[1,i] = ustrip(u"kpc", a[i].x)
            @inbounds Pos[2,i] = ustrip(u"kpc", a[i].y)
            @inbounds Pos[3,i] = ustrip(u"kpc", a[i].z)
        end
    elseif mode == "SI"
        for i in 1:NumTotal
            @inbounds Pos[1,i] = ustrip(u"m", a[i].x)
            @inbounds Pos[2,i] = ustrip(u"m", a[i].y)
            @inbounds Pos[3,i] = ustrip(u"m", a[i].z)
        end
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
    return kdtree = KDTree(Pos)
end

function kdtree_setup(a::Array{PhysicalParticle2D,1}; mode = "Astro")
    NumTotal = length(a)
    Pos = zeros(2,NumTotal)
    if mode == "Astro"
        for i in 1:NumTotal
            @inbounds Pos[1,i] = ustrip(u"kpc", a[i].Pos.x)
            @inbounds Pos[2,i] = ustrip(u"kpc", a[i].Pos.y)
        end
    elseif mode == "SI"
        for i in 1:NumTotal
            @inbounds Pos[1,i] = ustrip(u"m", a[i].Pos.x)
            @inbounds Pos[2,i] = ustrip(u"m", a[i].Pos.y)
        end
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
    return kdtree = KDTree(Pos)
end

function kdtree_setup(a::Array{PhysicalParticle3D,1}; mode = "Astro")
    NumTotal = length(a)
    Pos = zeros(3,NumTotal)
    if mode == "Astro"
        for i in 1:NumTotal
            @inbounds Pos[1,i] = ustrip(u"kpc", a[i].Pos.x)
            @inbounds Pos[2,i] = ustrip(u"kpc", a[i].Pos.y)
            @inbounds Pos[3,i] = ustrip(u"kpc", a[i].Pos.z)
        end
    elseif mode == "SI"
        for i in 1:NumTotal
            @inbounds Pos[1,i] = ustrip(u"m", a[i].Pos.x)
            @inbounds Pos[2,i] = ustrip(u"m", a[i].Pos.y)
            @inbounds Pos[3,i] = ustrip(u"m", a[i].Pos.z)
        end
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
    return kdtree = KDTree(Pos)
end

#####  K-nearest Search  #####

function kdtree_k_search(kdtree::KDTree, Point::Point2D, k::Int64; LeafSize = 10, SortByDistance = true)
    return IDs, Distances = knn(kdtree, [Point.x, Point.y], k, SortByDistance)
end

function kdtree_k_search(kdtree::KDTree, Point::Point3D, k::Int64; LeafSize = 10, SortByDistance = true)
    return IDs, Distances = knn(kdtree, [Point.x, Point.y, Point.z], k, SortByDistance)
end

function kdtree_k_search(kdtree::KDTree, Point::PhysicalVector2D, k::Int64; LeafSize = 10, SortByDistance = true, mode = "Astro")
    if mode == "Astro"
        return IDs, Distances = knn(kdtree, [ustrip(u"kpc", Point.x), ustrip(u"kpc", Point.y)], k, SortByDistance)
    elseif mode == "SI"
        return IDs, Distances = knn(kdtree, [ustrip(u"m", Point.x), ustrip(u"m", Point.y)], k, SortByDistance)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function kdtree_k_search(kdtree::KDTree, Point::PhysicalVector3D, k::Int64; LeafSize = 10, SortByDistance = true, mode = "Astro")
    if mode == "Astro"
        return IDs, Distances = knn(kdtree, [ustrip(u"kpc", Point.x), ustrip(u"kpc", Point.y), ustrip(u"kpc", Point.z)], k, SortByDistance)
    elseif mode == "SI"
        return IDs, Distances = knn(kdtree, [ustrip(u"m", Point.x), ustrip(u"m", Point.y), ustrip(u"m", Point.z)], k, SortByDistance)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function kdtree_k_search(kdtree::KDTree, x::Float64, y::Float64, k::Int64; LeafSize = 10, SortByDistance = true)
    return IDs, Distances = knn(kdtree, [x, y], k, SortByDistance)
end

function kdtree_k_search(kdtree::KDTree, x::Float64, y::Float64, z::Float64, k::Int64; LeafSize = 10, SortByDistance = true)
    return IDs, Distances = knn(kdtree, [x, y, z], k, SortByDistance)
end

function kdtree_k_search(kdtree::KDTree, x::Quantity, y::Quantity, k::Int64; LeafSize = 10, SortByDistance = true, mode = "Astro")
    if mode == "Astro"
        return IDs, Distances = knn(kdtree, [ustrip(u"kpc", x), ustrip(u"kpc", y)], k, SortByDistance)
    elseif mode == "SI"
        return IDs, Distances = knn(kdtree, [ustrip(u"m", x), ustrip(u"m", y)], k, SortByDistance)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function kdtree_k_search(kdtree::KDTree, x::Quantity, y::Quantity, z::Quantity, k::Int64; LeafSize = 10, SortByDistance = true, mode = "Astro")
    if mode == "Astro"
        return IDs, Distances = knn(kdtree, [ustrip(u"kpc", x), ustrip(u"kpc", y), ustrip(u"kpc", z)], k, SortByDistance)
    elseif mode == "SI"
        return IDs, Distances = knn(kdtree, [ustrip(u"m", x), ustrip(u"m", y), ustrip(u"m", z)], k, SortByDistance)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function kdtree_k_search(kdtree::KDTree, Point::Array{Float64,1}, k::Int64; LeafSize = 10, SortByDistance = true)
    return IDs, Distances = knn(kdtree, Point, k, SortByDistance)
end

function kdtree_k_search(kdtree::KDTree, Point::Array{Quantity,1}, k::Int64; LeafSize = 10, SortByDistance = true, mode = "Astro")
    if mode == "Astro"
        return IDs, Distances = knn(kdtree, ustrip.(u"kpc", p), k, SortByDistance)
    elseif mode == "SI"
        return IDs, Distances = knn(kdtree, ustrip.(u"m", p), k, SortByDistance)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

#####   Radius Search   #####

function kdtree_radius_search(kdtree::KDTree, Point::Point2D, Radius::Float64; LeafSize = 10, SortByDistance = true)
    return IDs, Distances = knn(kdtree, [Point.x, Point.y], Radius, SortByDistance)
end

function kdtree_radius_search(kdtree::KDTree, Point::Point3D, Radius::Float64; LeafSize = 10, SortByDistance = true)
    return IDs, Distances = knn(kdtree, [Point.x, Point.y, Point.z], Radius, SortByDistance)
end

function kdtree_radius_search(kdtree::KDTree, Point::PhysicalVector2D, Radius::Quantity; LeafSize = 10, SortByDistance = true, mode = "Astro")
    if mode == "Astro"
        return IDs, Distances = knn(kdtree, [ustrip(u"kpc", Point.x), ustrip(u"kpc", Point.y)], ustrip(u"kpc", Radius), SortByDistance)
    elseif mode == "SI"
        return IDs, Distances = knn(kdtree, [ustrip(u"m", Point.x), ustrip(u"m", Point.y)], ustrip(u"m", Radius), SortByDistance)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function kdtree_radius_search(kdtree::KDTree, Point::PhysicalVector3D, Radius::Quantity; LeafSize = 10, SortByDistance = true, mode = "Astro")
    if mode == "Astro"
        return IDs, Distances = knn(kdtree, [ustrip(u"kpc", Point.x), ustrip(u"kpc", Point.y), ustrip(u"kpc", Point.z)], ustrip(u"kpc", Radius), SortByDistance)
    elseif mode == "SI"
        return IDs, Distances = knn(kdtree, [ustrip(u"m", Point.x), ustrip(u"m", Point.y), ustrip(u"m", Point.z)], ustrip(u"m", Radius), SortByDistance)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function kdtree_radius_search(kdtree::KDTree, x::Float64, y::Float64, Radius::Float64; LeafSize = 10, SortByDistance = true)
    return IDs, Distances = knn(kdtree, [x, y], Radius, SortByDistance)
end

function kdtree_radius_search(kdtree::KDTree, x::Float64, y::Float64, z::Float64, Radius::Float64; LeafSize = 10, SortByDistance = true)
    return IDs, Distances = knn(kdtree, [x, y, z], Radius, SortByDistance)
end

function kdtree_radius_search(kdtree::KDTree, x::Quantity, y::Quantity, Radius::Quantity; LeafSize = 10, SortByDistance = true, mode = "Astro")
    if mode == "Astro"
        return IDs, Distances = knn(kdtree, [ustrip(u"kpc", x), ustrip(u"kpc", y)], ustrip(u"kpc", Radius), SortByDistance)
    elseif mode == "SI"
        return IDs, Distances = knn(kdtree, [ustrip(u"m", x), ustrip(u"m", y)], ustrip(u"m", Radius), SortByDistance)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function kdtree_radius_search(kdtree::KDTree, x::Quantity, y::Quantity, z::Quantity, Radius::Quantity; LeafSize = 10, SortByDistance = true, mode = "Astro")
    if mode == "Astro"
        return IDs, Distances = knn(kdtree, [ustrip(u"kpc", x), ustrip(u"kpc", y), ustrip(u"kpc", z)], ustrip(u"kpc", Radius), SortByDistance)
    elseif mode == "SI"
        return IDs, Distances = knn(kdtree, [ustrip(u"m", x), ustrip(u"m", y), ustrip(u"m", z)], ustrip(u"m", Radius), SortByDistance)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function kdtree_radius_search(kdtree::KDTree, Point::Array{Float64,1}, Radius::Float64; LeafSize = 10, SortByDistance = true)
    return IDs, Distances = knn(kdtree, Point, Radius, SortByDistance)
end

function kdtree_radius_search(kdtree::KDTree, Point::Array{Quantity,1}, Radius::Quantity; LeafSize = 10, SortByDistance = true, mode = "Astro")
    if mode == "Astro"
        return IDs, Distances = knn(kdtree, ustrip.(u"kpc", p), ustrip(u"kpc", Radius), SortByDistance)
    elseif mode == "SI"
        return IDs, Distances = knn(kdtree, ustrip.(u"m", p), ustrip(u"m", Radius), SortByDistance)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end
