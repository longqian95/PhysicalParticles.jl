include("../src/PhysicalParticles.jl")
using .PhysicalParticles
using Unitful, UnitfulAstro

a = rand(3)
p = PhysicalParticles.pconvert(a, u"m")
@show p

b = rand(3, 10)
bv = pconvert(b, u"m")

c = rand(3, 10)
cv = pconvert(c, u"m/s")
@show bv + cv*10.0u"s"
@show mean(bv)
