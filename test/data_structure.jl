# include("./src/PhysicalParticles.jl"); using .PhysicalParticles, Unitful, UnitfulAstro

NumProcessors = 4
NumTotal = 5000
NumGas = 5000

using Distributed
using BenchmarkTools

addprocs(NumProcessors)
@everywhere using DistributedArrays, Unitful, UnitfulAstro
@everywhere using PhysicalParticles

@everywhere mutable struct ParticleData
    particles::Array{PhysicalParticle3D,1}
    gases::Array{GasParticle3D,1}
end

@everywhere mutable struct UniversalParticle
    Pos::PhysicalVector3D
    Vel::PhysicalVector3D
    Acc::PhysicalVector3D
    Mass::Quantity
    ID::Int64
    Type::ParticleType

    Entropy::Quantity
    Density::Quantity
    Hsml::Quantity

    RotVel::PhysicalVector3D
    DivVel::Quantity
    CurlVel::Quantity
    dHsmlRho::Float64

    Pressure::Quantity
    DtEntropy::Quantity
    MaxSignalVel::Quantity
end
@everywhere UniversalParticle() = UniversalParticle(PositionAstro(0.0,0.0,0.0), VelocityAstro(0.0,0.0,0.0), AccelerationAstro(0.0,0.0,0.0),
                                        0.0u"Msun", 0, ParticleType(5),
                                        0.0u"J/K", 0.0u"Msun/kpc^3", 0.0u"kpc",
                                        VelocityAstro(0.0,0.0,0.0), 0.0u"Gyr^-1", 0.0u"Gyr^-1", 0.0,
                                        0.0u"N/m^2", 0.0u"J/K/s", 0.0u"kpc/s")

@everywhere function move(a::DArray)
    for p in a
        p.Pos += p.Vel * 1.0u"Gyr"
    end
end

@everywhere function move(a::UniversalParticle)
    for p in a
        if p.Type == ParticleType(5)
            p.Pos += p.Vel * 1.0u"Gyr"
        elseif p.Type == ParticleType(1)
            p.Pos += p.Vel * 1.0u"Gyr"
        end
    end
end

@everywhere function move(a::Array{AbstractParticle3D})
    for p in a
        if typeof(p) == PhysicalParticle3D
            p.Pos += p.Vel * 1.0u"Gyr"
        elseif typeof(p) == GasParticle3D
            p.Pos += p.Vel * 1.0u"Gyr"
        end

    end
end

##### 1 Different Arrays for different particle types #####

function f1()
    array1 = ParticleData(fill(PhysicalParticle3D(),NumTotal), fill(GasParticle3D(),NumGas))
    Darray1_particle = distribute(array1.particles)
    Darray1_gas = distribute(array1.gases)
    move(Darray1_particle)
    move(Darray1_gas)
end
t1 = @benchmarkable f1()
test1() = run(t1, samples = 10)

##### 2 Universal particle type #####

function f2()
    array2 = fill(UniversalParticle(), NumTotal)
    Darray2 = distribute(array2)
    move(Darray2)
end
t2 = @benchmarkable f2()
test2() = run(t2, samples = 10)

##### 3 Abstract type #####

function f3()
    array3 = [fill(PhysicalParticle3D(),NumTotal); fill(GasParticle3D(),NumGas)]
    Darray3 = distribute(array3)
    move(Darray3)
end
t3 = @benchmarkable f3()
test3() = run(t3, samples = 10)

#println("#####  Test 1 : Different Arrays for different particle types  #####")
#test1()
#
#println("#####  Test 2 : Universal particle type  #####")
#test2()
#
#println("#####  Test 3 : Abstract type  #####")
#test3()
