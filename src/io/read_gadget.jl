"""
    function read_gadget2(filename::String)
        Read Gadget2 format snapshots and return a struct of StarData
        Returns a turple:
            (::Header_Gadget2, ::Array{PhysicalParticle,1}, ::Array{GasData,1})
"""
function read_gadget2(filename::String)
    f = open(filename, "r")

    # Read Header
    temp1 = read(f, Int32)
    Header = Header_Gadget2()

    for i in 1:6
        Header.npart[i] = read(f, Int32)
    end
    for i in 1:6
        Header.mass[i] = read(f, Float64)
    end
    Header.time = read(f, Float64)
    Header.redshift = read(f, Float64)
    Header.flag_sfr = read(f, Int32)
    Header.flag_feedback = read(f, Int32)
    for i in 1:6
        Header.npartTotal[i] = read(f, UInt32)
    end
    Header.flag_cooling = read(f, Int32)
    Header.num_files = read(f, Int32)
    Header.BoxSize = read(f, Float64)
    Header.Omega0 = read(f, Float64)
    Header.OmegaLambda = read(f, Float64)
    Header.HubbleParam = read(f, Float64)
    Header.flag_stellarage = read(f, Int32)
    Header.flag_metals = read(f, Int32)
    for i in 1:6
        Header.npartTotalHighWord[i] = read(f, UInt32)
    end
    Header.flag_entropy_instead_u = read(f, Int32)

    NumTotal = sum(Header.npart)
    Particles = fill(PhysicalParticle(), NumTotal)

    for i = 1:60
        read(f, 1)
    end
    temp2 = read(f, Int32)
    if temp1!=temp2
        error("Wrong location symbol!\n")
        quit()
    end

    # Read Position Block
    temp1 = read(f, Int32)
    for i in 1:NumTotal
        x = read(f, Float32)*u"kpc"
        y = read(f, Float32)*u"kpc"
        z = read(f, Float32)*u"kpc"
        Particles[i].Pos = PhysicalVector(x,y,z)
    end
    temp2 = read(f, Int32)
    if temp1!=temp2
        error("Wrong location symbol!\n")
        quit()
    end

    # Read Velocity Block
    temp1 = read(f, Int32)
    Vel = zeros(3,NumTotal)
    for i in 1:NumTotal
        vx = read(f, Float32)*u"kpc/Gyr"
        vy = read(f, Float32)*u"kpc/Gyr"
        vz = read(f, Float32)*u"kpc/Gyr"
        Particles[i].Vel = PhysicalVector(vx,vy,vz)
    end
    temp2 = read(f, Int32)
    if temp1!=temp2
        error("Wrong location symbol!\n")
        quit()
    end

    # Read ID Block
    temp1 = read(f, Int32)
    ID = zeros(UInt32, NumTotal)
    for i in 1:NumTotal
        Particles[i].ID = read(f, UInt32)
    end
    temp2 = read(f, Int32)
    if temp1!=temp2
        error("Wrong location symbol!\n")
        quit()
    end

    # Read Mass Block
    if sum(Header.mass)==0.0
        temp1 = read(f, Int32)
    end

    count_temp = 0
    for i in 1:6
        for k in (count_temp+1):(count_temp+Header.npart[i])
            if Header.mass[i] == 0.0
                Particles[i].Mass = read(f, Float32)*u"Msun"
            else
                Particles[i].Mass = Header.mass[i]*u"Msun"
            end
        end
        count_temp += Header.npart[i]
    end

    if sum(Header.mass)==0.0
        temp2 = read(f, Int32)
        if temp1!=temp2
            error("Wrong location symbol!\n")
            quit()
        end
    end

    # Read Gas Internal Energy Block
    NumGas = Header.npart[1]
    SphData = fill(GasData(), NumGas)
    if NumGas>0
        # Read Entropy
        temp1 = read(f, Int32)
        for i in 1:NumGas
            SphData[i].Entropy = read(f, Float32)*u"J/K"
        end
        temp2 = read(f, Int32)
        if temp1!=temp2
            error("Wrong location symbol!\n")
            quit()
        end

        # Read Density
        if !eof(f)
            temp1 = read(f, Int32)
            for i in 1:NumGas
                SphData[i].Density = read(f, Float32)*u"Msun/kpc^3"
            end
            temp2 = read(f, Int32)
            if temp1!=temp2
                error("Wrong location symbol!\n")
                quit()
            end
        end

        # Read Hsml
        if !eof(f)
            temp1 = read(f, Int32)
            for i in 1:NumGas
                SphData[i].Hsml = read(f, Float32)*u"kpc"
            end
            temp2 = read(f, Int32)
            if temp1!=temp2
                error("Wrong location symbol!\n")
                quit()
            end
        end
    end

    return Header, Particles, SphData
end
