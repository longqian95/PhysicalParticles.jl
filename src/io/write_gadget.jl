"""
    function write_gadget2(filename::String, ParticleData::StarData)
        Saves data in Gadget2 format
        Header.npart[1] determines whether to output gasous data blocks
"""
function write_gadget2(filename::String, Header::Header_Gadget2, Particles::Array{PhysicalParticle,1}, SphData::Array{GasData,1})
    f = open(filename, "w")

    # Write Header
    temp::Int32 = 256
    write(f, temp)
    for i in 1:6
        write(f, Int32(Header.npart[i]))
    end
    for i in 1:6
        write(f, Float64(Header.mass[i]))
    end
    write(f, Float64(Header.time)           )
    write(f, Float64(Header.redshift)       )
    write(f, Int32(Header.flag_sfr)       )
    write(f, Int32(Header.flag_feedback)  )
    for i in 1:6
        write(f, UInt32(Header.npartTotal[i]))
    end
    write(f, Int32(Header.flag_cooling)   )
    write(f, Int32(Header.num_files)      )
    write(f, Float64(Header.BoxSize)        )
    write(f, Float64(Header.Omega0)         )
    write(f, Float64(Header.OmegaLambda)    )
    write(f, Float64(Header.HubbleParam)    )
    write(f, Int32(Header.flag_stellarage))
    write(f, Int32(Header.flag_metals)    )
    for i in 1:6
        write(f, UInt32(Header.npartTotalHighWord[i]))
    end
    write(f, Int32(Header.flag_entropy_instead_u))
    for i = 1:60
        write(f, "\0")
    end
    write(f, temp)

    NumTotal = length(Particles)
    NumGas = length(SphData)

    # Write Position Block
    temp = 4 * NumTotal * 3
    write(f, temp)
    for i in 1:NumTotal
        write(f, Float32(ustrip(u"kpc", Particles[i].Pos.x)))
        write(f, Float32(ustrip(u"kpc", Particles[i].Pos.y)))
        write(f, Float32(ustrip(u"kpc", Particles[i].Pos.z)))
    end
    write(f, temp)

    # Write Velocity Block
    write(f, temp)
    for i in 1:NumTotal
        write(f, Float32(ustrip(u"kpc/Gyr", Particles[i].Vel.x)))
        write(f, Float32(ustrip(u"kpc/Gyr", Particles[i].Vel.y)))
        write(f, Float32(ustrip(u"kpc/Gyr", Particles[i].Vel.z)))
    end
    write(f, temp)

    # Write ID Block
    temp = 4 * NumTotal
    write(f, temp)
    for i in 1:NumTotal
        write(f, Particles[i].ID)
    end
    write(f, temp)

    # Write Mass Block
    if sum(Header.mass)==0.0
        write(f, temp)
    end

    count_temp = 0
    for i in 1:6
        for k in (count_temp+1):(count_temp+Header.npart[i])
            if Header.mass[i] == 0.0
                write(f, Float32(ustrip(u"Msun", Particles[i].Mass)))
            end
        end
        count_temp += Header.npart[i]
    end

    if sum(Header.mass)==0.0
        write(f, temp)
    end

    # Write Gas
    NumGas = Header.npart[1]
    if NumGas>0
        # Write Entropy
        write(f, temp)
        for i in 1:NumGas
            write(f, Float32(ustrip(u"J/K", SphData[i].Entropy)))
        end
        write(f, temp)

        # Write Density
        write(f, temp)
        for i in 1:NumGas
            write(f, Float32(ustrip(u"Msun/kpc^3", SphData[i].Density)))
        end
        write(f, temp)

        # Write Hsml
        write(f, temp)
        for i in 1:NumGas
            write(f, Float32(ustrip(u"kpc", SphData[i].Hsml)))
        end
        write(f, temp)
    end

    close(f)
    return true
end
