function write_ascii(FileName::String, array::Array{Point2D})
    file = open(FileName, "w")
    write(file, "#id x y\n")
    for i in 1:length(array)
        @inbounds write(file, "$i $(array[i].x) $(array[i].y)\n")
    end
    close(file)
end

function write_ascii(FileName::String, array::Array{Point3D})
    file = open(FileName, "w")
    write(file, "#id x y z\n")
    for i in 1:length(array)
        @inbounds write(file, "$i $(array[i].x) $(array[i].y) $(array[i].z))\n")
    end
    close(file)
end

function write_ascii(FileName::String, array::Array{PhysicalVector2D}; mode = "Astro")
    file = open(FileName, "w")
    if mode == "Astro"
        write(file, "#id x/[kpc] y/[kpc]\n")
        for i in 1:length(array)
            @inbounds write(file, "$i $(ustrip(UnitfulAstro.kpc, array[i].x)) $(ustrip(UnitfulAstro.kpc, array[i].y))\n")
        end
        close(file)
    elseif mode == "SI"
        write(file, "#id x/[m] y/[m]\n")
        for i in 1:length(array)
            @inbounds write(file, "$i $(ustrip(UnitfulAstro.m, array[i].x)) $(ustrip(UnitfulAstro.m, array[i].y))\n")
        end
        close(file)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function write_ascii(FileName::String, array::Array{PhysicalVector3D}; mode = "Astro")
    file = open(FileName, "w")
    if mode == "Astro"
        write(file, "#id x/[kpc] z/[kpc]\n")
        for i in 1:length(array)
            @inbounds write(file, "$i $(ustrip(UnitfulAstro.kpc, array[i].x)) $(ustrip(UnitfulAstro.kpc, array[i].y)) $(ustrip(UnitfulAstro.kpc, array[i].z))\n")
        end
        close(file)
    elseif mode == "SI"
        write(file, "#id x/[m] y/[m] z/[m]\n")
        for i in 1:length(array)
            @inbounds write(file, "$i $(ustrip(UnitfulAstro.m, array[i].x)) $(ustrip(UnitfulAstro.m, array[i].y)) $(ustrip(UnitfulAstro.m, array[i].z))\n")
        end
        close(file)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function write_ascii(FileName::String, array::Array{PhysicalParticle2D}; mode = "Astro")
    file = open(FileName, "w")
    if mode == "Astro"
        write(file, "#id x y [kpc] || vx vy [kpc/Gyr] || ax ay [kpc/Gyr^2] || m [Msun] || Type\n")
        for i in 1:length(array)
            @inbounds write(file, "$i $(ustrip(UnitfulAstro.kpc, array[i].Pos.x)) $(ustrip(UnitfulAstro.kpc, array[i].Pos.y)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.kpc/UnitfulAstro.Gyr, array[i].Vel.x)) $(ustrip(UnitfulAstro.kpc/UnitfulAstro.Gyr, array[i].Vel.y)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.kpc/UnitfulAstro.Gyr^2, array[i].Acc.x)) $(ustrip(UnitfulAstro.kpc/UnitfulAstro.Gyr^2, array[i].Acc.y)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.Msun, array[i].Mass)) $(array[i].Type)\n")
        end
        close(file)
    elseif mode == "SI"
        write(file, "#id x y [m] || vx vy [m/s] || ax ay [m/s^2] || m [kg] || Type\n")
        for i in 1:length(array)
            @inbounds write(file, "$i $(ustrip(UnitfulAstro.m, array[i].Pos.x)) $(ustrip(UnitfulAstro.m, array[i].Pos.y)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.m/UnitfulAstro.s, array[i].Vel.x)) $(ustrip(UnitfulAstro.m/UnitfulAstro.s, array[i].Vel.y)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.m/UnitfulAstro.s^2, array[i].Acc.x)) $(ustrip(UnitfulAstro.m/UnitfulAstro.s^2, array[i].Acc.y)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.kg, array[i].Mass)) $(array[i].Type)\n")
        end
        close(file)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end

function write_ascii(FileName::String, array::Array{PhysicalParticle3D}; mode = "Astro")
    file = open(FileName, "w")
    if mode == "Astro"
        write(file, "#id x y z [kpc] || vx vy vz [kpc/Gyr] || ax ay az [kpc/Gyr^2] || m [Msun] || Type\n")
        for i in 1:length(array)
            @inbounds write(file, "$i $(ustrip(UnitfulAstro.kpc, array[i].Pos.x)) $(ustrip(UnitfulAstro.kpc, array[i].Pos.y)) $(ustrip(UnitfulAstro.kpc, array[i].Pos.z)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.kpc/UnitfulAstro.Gyr, array[i].Vel.x)) $(ustrip(UnitfulAstro.kpc/UnitfulAstro.Gyr, array[i].Vel.y)) $(ustrip(UnitfulAstro.kpc/UnitfulAstro.Gyr, array[i].Vel.z)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.kpc/UnitfulAstro.Gyr^2, array[i].Acc.x)) $(ustrip(UnitfulAstro.kpc/UnitfulAstro.Gyr^2, array[i].Acc.y)) $(ustrip(UnitfulAstro.kpc/UnitfulAstro.Gyr^2, array[i].Acc.z)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.Msun, array[i].Mass)) $(array[i].Type)\n")
        end
        close(file)
    elseif mode == "SI"
        write(file, "#id x y z [m] || vx vy vz [m/s] || ax ay az [m/s^2] || m [kg] || Type\n")
        for i in 1:length(array)
            @inbounds write(file, "$i $(ustrip(UnitfulAstro.m, array[i].Pos.x)) $(ustrip(UnitfulAstro.m, array[i].Pos.y)) $(ustrip(UnitfulAstro.m, array[i].Pos.z)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.m/UnitfulAstro.s, array[i].Vel.x)) $(ustrip(UnitfulAstro.m/UnitfulAstro.s, array[i].Vel.y)) $(ustrip(UnitfulAstro.m/UnitfulAstro.s, array[i].Vel.z)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.m/UnitfulAstro.s^2, array[i].Acc.x)) $(ustrip(UnitfulAstro.m/UnitfulAstro.s^2, array[i].Acc.y)) $(ustrip(UnitfulAstro.m/UnitfulAstro.s^2, array[i].Acc.z)) ")
            @inbounds write(file, "$(ustrip(UnitfulAstro.kg, array[i].Mass)) $(array[i].Type)\n")
        end
        close(file)
    else
        error("Not supported unit mode!\n    Supported: Astro (default) and SI")
    end
end
