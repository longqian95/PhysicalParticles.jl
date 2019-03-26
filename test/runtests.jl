include("../src/PhysicalParticles.jl")

using .PhysicalParticles
using Test, Unitful, UnitfulAstro

array2d = [1.0 2.0 3.0 4.0 5.0;
           6.0 7.0 8.0 9.0 10.0]
array3d = [1.0 2.0 3.0 4.0 5.0;
           6.0 7.0 8.0 9.0 10.0;
           11.0 12.0 13.0 14.0 15.0]


@testset "2D Points" begin
    a = Point(1.0,2.0)
    b = Point(2.0,2.0)
    @test norm(a+b) == 5.0
    @test norm(b-a) == 1.0
    @test a*b == 6.0
    @test dot(a,b) == 6.0
    @test (a+1.0)*(1.0+b) == 15.0
    @test (a-1.0)*(1.0-b) == -1.0
    @test gety(a*2.0) == 4.0
    @test getx(a/2.0) == 0.5
    @test (2.0*a).x == 2.0
    @test zero(a).y == 0.0

    @test normalize(Point(3.0,4.0)).x == 0.6
    @test rotate(Point(1.0,0.0), 0.5pi).y == 1.0
    @test norm(npconvert([3.0,4.0])) == 5.0

    c = npconvert(array2d)
    @test mean(c) == Point(3.0,8.0)
    @test center(c) == Point(3.0,8.0)
    @test min_x(c) == 1.0
    @test max_x(c) == 5.0
    @test min_y(c) == 6.0
    @test max_y(c) == 10.0
end

@testset "3D Points" begin
    a = Point(1.0,2.0,3.0)
    b = Point(2.0,2.0,3.0)
    @test a+b == Point(3.0,4.0,6.0)
    @test norm(b-a) == 1.0
    @test a*b == 15.0
    @test dot(a,b) == 15.0
    @test (a+1.0)*(1.0+b) == 31.0
    @test (a-1.0)*(1.0-b) == -5.0
    @test gety(a*2.0) == 4.0
    @test getx(a/2.0) == 0.5
    @test getz(2.0*a) == 6.0
    @test zero(a).y == 0.0

    @test normalize(Point(3.0,4.0,12.0)).x == 3.0/13.0
    @test rotate_x(Point(0.0,1.0,0.0), 0.5pi).z == 1.0
    @test rotate_y(Point(0.0,0.0,1.0), 0.5pi).x == 1.0
    @test rotate_z(Point(1.0,0.0,0.0), 0.5pi).y == 1.0
    @test norm(npconvert([3.0,4.0,12.0])) == 13.0

    c = npconvert(array3d)
    @test mean(c) == Point(3.0,8.0,13.0)
    @test center(c) == Point(3.0,8.0,13.0)
    @test min_x(c) == 1.0
    @test max_x(c) == 5.0
    @test min_y(c) == 6.0
    @test max_y(c) == 10.0
    @test min_z(c) == 11.0
    @test max_z(c) == 15.0
end

@testset "2D Physical Vector" begin
    a = PhysicalVector(1.0,2.0,u"m")
    b = PhysicalVector(2.0,2.0,u"m")
    @test norm(a+b) == 5.0u"m"
    @test norm(b-a) == 1.0u"m"
    @test a*b == 6.0u"m^2"
    @test dot(a,b) == 6.0u"m^2"
    @test (a+1.0u"m")*(1.0u"m"+b) == 15.0u"m^2"
    @test (a-1.0u"m")*(1.0u"m"-b) == -1.0u"m^2"
    @test gety(a*2.0) == 4.0u"m"
    @test getx(a/2.0) == 0.5u"m"
    @test (2.0*a).x == 2.0u"m"
    @test zero(a).y == 0.0u"m"

    @test normalize(PhysicalVector(3.0,4.0,u"m")).x == 0.6u"m"
    @test rotate(PhysicalVector(1.0,0.0,u"m"), 0.5pi).y == 1.0u"m"
    @test norm(pconvert([3.0,4.0],u"m")) == 5.0u"m"

    c = pconvert(array2d,u"m")
    @test mean(c) == PhysicalVector(3.0,8.0,u"m")
    @test center(c) == PhysicalVector(3.0,8.0,u"m")
    @test min_x(c) == 1.0u"m"
    @test max_x(c) == 5.0u"m"
    @test min_y(c) == 6.0u"m"
    @test max_y(c) == 10.0u"m"
end

@testset "2D Position" begin
    a = Position(1.0,2.0)
    b = Position(2.0,2.0)
    @test norm(a+b) == 5.0u"m"
    @test norm(b-a) == 1.0u"m"
    @test a*b == 6.0u"m^2"
    @test dot(a,b) == 6.0u"m^2"
    @test (a+1.0u"m")*(1.0u"m"+b) == 15.0u"m^2"
    @test (a-1.0u"m")*(1.0u"m"-b) == -1.0u"m^2"
    @test gety(a*2.0) == 4.0u"m"
    @test getx(a/2.0) == 0.5u"m"
    @test (2.0*a).x == 2.0u"m"
    @test zero(a).y == 0.0u"m"

    @test normalize(Position(3.0,4.0)).x == 0.6u"m"
    @test rotate(Position(1.0,0.0), 0.5pi).y == 1.0u"m"
    @test norm(pconvert([3.0,4.0],u"m")) == 5.0u"m"

    c = pconvert(array2d,u"m")
    @test mean(c) == Position(3.0,8.0)
    @test center(c) == Position(3.0,8.0)
    @test min_x(c) == 1.0u"m"
    @test max_x(c) == 5.0u"m"
    @test min_y(c) == 6.0u"m"
    @test max_y(c) == 10.0u"m"
end

@testset "2D Position Astro" begin
    a = PositionAstro(1.0,2.0)
    b = PositionAstro(2.0,2.0)
    @test norm(a+b) == 5.0u"kpc"
    @test norm(b-a) == 1.0u"kpc"
    @test a*b == 6.0u"kpc^2"
    @test dot(a,b) == 6.0u"kpc^2"
    @test (a+1.0u"kpc")*(1.0u"kpc"+b) == 15.0u"kpc^2"
    @test (a-1.0u"kpc")*(1.0u"kpc"-b) == -1.0u"kpc^2"
    @test gety(a*2.0) == 4.0u"kpc"
    @test getx(a/2.0) == 0.5u"kpc"
    @test (2.0*a).x == 2.0u"kpc"
    @test zero(a).y == 0.0u"kpc"

    @test normalize(PositionAstro(3.0,4.0)).x == 0.6u"kpc"
    @test rotate(PositionAstro(1.0,0.0), 0.5pi).y == 1.0u"kpc"
    @test norm(pconvert([3.0,4.0],u"kpc")) == 5.0u"kpc"

    c = pconvert(array2d,u"kpc")
    @test mean(c) == PositionAstro(3.0,8.0)
    @test center(c) == PositionAstro(3.0,8.0)
    @test min_x(c) == 1.0u"kpc"
    @test max_x(c) == 5.0u"kpc"
    @test min_y(c) == 6.0u"kpc"
    @test max_y(c) == 10.0u"kpc"
end

@testset "2D Velocity" begin
    a = Velocity(1.0,2.0)
    b = Velocity(2.0,2.0)
    @test norm(a+b) == 5.0u"m/s"
    @test norm(b-a) == 1.0u"m/s"
    @test a*b == 6.0u"m^2/s^2"
    @test dot(a,b) == 6.0u"m^2/s^2"
    @test (a+1.0u"m/s")*(1.0u"m/s"+b) == 15.0u"m^2/s^2"
    @test (a-1.0u"m/s")*(1.0u"m/s"-b) == -1.0u"m^2/s^2"
    @test gety(a*2.0) == 4.0u"m/s"
    @test getx(a/2.0) == 0.5u"m/s"
    @test (2.0*a).x == 2.0u"m/s"
    @test zero(a).y == 0.0u"m/s"

    @test normalize(Velocity(3.0,4.0)).x == 0.6u"m/s"
    @test rotate(Velocity(1.0,0.0), 0.5pi).y == 1.0u"m/s"
    @test norm(pconvert([3.0,4.0],u"m/s")) == 5.0u"m/s"

    c = pconvert(array2d,u"m/s")
    @test mean(c) == Velocity(3.0,8.0)
    @test center(c) == Velocity(3.0,8.0)
    @test min_x(c) == 1.0u"m/s"
    @test max_x(c) == 5.0u"m/s"
    @test min_y(c) == 6.0u"m/s"
    @test max_y(c) == 10.0u"m/s"
end

@testset "2D Velocity Astro" begin
    a = VelocityAstro(1.0,2.0)
    b = VelocityAstro(2.0,2.0)
    @test norm(a+b) == 5.0u"kpc/Gyr"
    @test norm(b-a) == 1.0u"kpc/Gyr"
    @test a*b == 6.0u"kpc^2/Gyr^2"
    @test dot(a,b) == 6.0u"kpc^2/Gyr^2"
    @test (a+1.0u"kpc/Gyr")*(1.0u"kpc/Gyr"+b) == 15.0u"kpc^2/Gyr^2"
    @test (a-1.0u"kpc/Gyr")*(1.0u"kpc/Gyr"-b) == -1.0u"kpc^2/Gyr^2"
    @test gety(a*2.0) == 4.0u"kpc/Gyr"
    @test getx(a/2.0) == 0.5u"kpc/Gyr"
    @test (2.0*a).x == 2.0u"kpc/Gyr"
    @test zero(a).y == 0.0u"kpc/Gyr"

    @test normalize(VelocityAstro(3.0,4.0)).x == 0.6u"kpc/Gyr"
    @test rotate(VelocityAstro(1.0,0.0), 0.5pi).y == 1.0u"kpc/Gyr"
    @test norm(pconvert([3.0,4.0],u"kpc/Gyr")) == 5.0u"kpc/Gyr"

    c = pconvert(array2d,u"kpc/Gyr")
    @test mean(c) == VelocityAstro(3.0,8.0)
    @test center(c) == VelocityAstro(3.0,8.0)
    @test min_x(c) == 1.0u"kpc/Gyr"
    @test max_x(c) == 5.0u"kpc/Gyr"
    @test min_y(c) == 6.0u"kpc/Gyr"
    @test max_y(c) == 10.0u"kpc/Gyr"
end

@testset "2D Acceleration" begin
    a = Acceleration(1.0,2.0)
    b = Acceleration(2.0,2.0)
    @test norm(a+b) == 5.0u"m/s^2"
    @test norm(b-a) == 1.0u"m/s^2"
    @test a*b == 6.0u"m^2/s^4"
    @test dot(a,b) == 6.0u"m^2/s^4"
    @test (a+1.0u"m/s^2")*(1.0u"m/s^2"+b) == 15.0u"m^2/s^4"
    @test (a-1.0u"m/s^2")*(1.0u"m/s^2"-b) == -1.0u"m^2/s^4"
    @test gety(a*2.0) == 4.0u"m/s^2"
    @test getx(a/2.0) == 0.5u"m/s^2"
    @test (2.0*a).x == 2.0u"m/s^2"
    @test zero(a).y == 0.0u"m/s^2"

    @test normalize(Acceleration(3.0,4.0)).x == 0.6u"m/s^2"
    @test rotate(Acceleration(1.0,0.0), 0.5pi).y == 1.0u"m/s^2"
    @test norm(pconvert([3.0,4.0],u"m/s^2")) == 5.0u"m/s^2"

    c = pconvert(array2d,u"m/s^2")
    @test mean(c) == Acceleration(3.0,8.0)
    @test center(c) == Acceleration(3.0,8.0)
    @test min_x(c) == 1.0u"m/s^2"
    @test max_x(c) == 5.0u"m/s^2"
    @test min_y(c) == 6.0u"m/s^2"
    @test max_y(c) == 10.0u"m/s^2"
end

@testset "2D Acceleration Astro" begin
    a = AccelerationAstro(1.0,2.0)
    b = AccelerationAstro(2.0,2.0)
    @test norm(a+b) == 5.0u"kpc/Gyr^2"
    @test norm(b-a) == 1.0u"kpc/Gyr^2"
    @test a*b == 6.0u"kpc^2/Gyr^4"
    @test dot(a,b) == 6.0u"kpc^2/Gyr^4"
    @test (a+1.0u"kpc/Gyr^2")*(1.0u"kpc/Gyr^2"+b) == 15.0u"kpc^2/Gyr^4"
    @test (a-1.0u"kpc/Gyr^2")*(1.0u"kpc/Gyr^2"-b) == -1.0u"kpc^2/Gyr^4"
    @test gety(a*2.0) == 4.0u"kpc/Gyr^2"
    @test getx(a/2.0) == 0.5u"kpc/Gyr^2"
    @test (2.0*a).x == 2.0u"kpc/Gyr^2"
    @test zero(a).y == 0.0u"kpc/Gyr^2"

    @test normalize(AccelerationAstro(3.0,4.0)).x == 0.6u"kpc/Gyr^2"
    @test rotate(AccelerationAstro(1.0,0.0), 0.5pi).y == 1.0u"kpc/Gyr^2"
    @test norm(pconvert([3.0,4.0],u"kpc/Gyr^2")) == 5.0u"kpc/Gyr^2"

    c = pconvert(array2d,u"kpc/Gyr^2")
    @test mean(c) == AccelerationAstro(3.0,8.0)
    @test center(c) == AccelerationAstro(3.0,8.0)
    @test min_x(c) == 1.0u"kpc/Gyr^2"
    @test max_x(c) == 5.0u"kpc/Gyr^2"
    @test min_y(c) == 6.0u"kpc/Gyr^2"
    @test max_y(c) == 10.0u"kpc/Gyr^2"
end

@testset "3D Physical Vector" begin
    a = PhysicalVector(1.0,2.0,3.0,u"m")
    b = PhysicalVector(2.0,2.0,3.0,u"m")
    @test a+b == PhysicalVector(3.0,4.0,6.0,u"m")
    @test norm(b-a) == 1.0u"m"
    @test a*b == 15.0u"m^2"
    @test dot(a,b) == 15.0u"m^2"
    @test (a+1.0u"m")*(1.0u"m"+b) == 31.0u"m^2"
    @test (a-1.0u"m")*(1.0u"m"-b) == -5.0u"m^2"
    @test gety(a*2.0) == 4.0u"m"
    @test getx(a/2.0) == 0.5u"m"
    @test getz(2.0*a) == 6.0u"m"
    @test zero(a).y == 0.0u"m"

    @test normalize(PhysicalVector(3.0,4.0,12.0,u"m")).x == 3.0/13.0*u"m"
    @test rotate_x(PhysicalVector(0.0,1.0,0.0,u"m"), 0.5pi).z == 1.0u"m"
    @test rotate_y(PhysicalVector(0.0,0.0,1.0,u"m"), 0.5pi).x == 1.0u"m"
    @test rotate_z(PhysicalVector(1.0,0.0,0.0,u"m"), 0.5pi).y == 1.0u"m"
    @test norm(pconvert([3.0,4.0,12.0],u"m")) == 13.0u"m"

    c = pconvert(array3d,u"m")
    @test mean(c) == PhysicalVector(3.0,8.0,13.0,u"m")
    @test center(c) == PhysicalVector(3.0,8.0,13.0,u"m")
    @test min_x(c) == 1.0u"m"
    @test max_x(c) == 5.0u"m"
    @test min_y(c) == 6.0u"m"
    @test max_y(c) == 10.0u"m"
    @test min_z(c) == 11.0u"m"
    @test max_z(c) == 15.0u"m"
end

@testset "3D Position" begin
    a = Position(1.0,2.0,3.0)
    b = Position(2.0,2.0,3.0)
    @test a+b == Position(3.0,4.0,6.0)
    @test norm(b-a) == 1.0u"m"
    @test a*b == 15.0u"m^2"
    @test dot(a,b) == 15.0u"m^2"
    @test (a+1.0u"m")*(1.0u"m"+b) == 31.0u"m^2"
    @test (a-1.0u"m")*(1.0u"m"-b) == -5.0u"m^2"
    @test gety(a*2.0) == 4.0u"m"
    @test getx(a/2.0) == 0.5u"m"
    @test getz(2.0*a) == 6.0u"m"
    @test zero(a).y == 0.0u"m"

    @test normalize(Position(3.0,4.0,12.0)).x == 3.0/13.0*u"m"
    @test rotate_x(Position(0.0,1.0,0.0), 0.5pi).z == 1.0u"m"
    @test rotate_y(Position(0.0,0.0,1.0), 0.5pi).x == 1.0u"m"
    @test rotate_z(Position(1.0,0.0,0.0), 0.5pi).y == 1.0u"m"
    @test norm(pconvert([3.0,4.0,12.0],u"m")) == 13.0u"m"

    c = pconvert(array3d,u"m")
    @test mean(c) == Position(3.0,8.0,13.0)
    @test center(c) == Position(3.0,8.0,13.0)
    @test min_x(c) == 1.0u"m"
    @test max_x(c) == 5.0u"m"
    @test min_y(c) == 6.0u"m"
    @test max_y(c) == 10.0u"m"
    @test min_z(c) == 11.0u"m"
    @test max_z(c) == 15.0u"m"
end

@testset "3D Position Astro" begin
    a = PositionAstro(1.0,2.0,3.0)
    b = PositionAstro(2.0,2.0,3.0)
    @test a+b == PositionAstro(3.0,4.0,6.0)
    @test norm(b-a) == 1.0u"kpc"
    @test a*b == 15.0u"kpc^2"
    @test dot(a,b) == 15.0u"kpc^2"
    @test (a+1.0u"kpc")*(1.0u"kpc"+b) == 31.0u"kpc^2"
    @test (a-1.0u"kpc")*(1.0u"kpc"-b) == -5.0u"kpc^2"
    @test gety(a*2.0) == 4.0u"kpc"
    @test getx(a/2.0) == 0.5u"kpc"
    @test getz(2.0*a) == 6.0u"kpc"
    @test zero(a).y == 0.0u"kpc"

    @test normalize(PositionAstro(3.0,4.0,12.0)).x == 3.0/13.0*u"kpc"
    @test rotate_x(PositionAstro(0.0,1.0,0.0), 0.5pi).z == 1.0u"kpc"
    @test rotate_y(PositionAstro(0.0,0.0,1.0), 0.5pi).x == 1.0u"kpc"
    @test rotate_z(PositionAstro(1.0,0.0,0.0), 0.5pi).y == 1.0u"kpc"
    @test norm(pconvert([3.0,4.0,12.0],u"kpc")) == 13.0u"kpc"

    c = pconvert(array3d,u"kpc")
    @test mean(c) == PositionAstro(3.0,8.0,13.0)
    @test center(c) == PositionAstro(3.0,8.0,13.0)
    @test min_x(c) == 1.0u"kpc"
    @test max_x(c) == 5.0u"kpc"
    @test min_y(c) == 6.0u"kpc"
    @test max_y(c) == 10.0u"kpc"
    @test min_z(c) == 11.0u"kpc"
    @test max_z(c) == 15.0u"kpc"
end

@testset "3D Velocity" begin
    a = Velocity(1.0,2.0,3.0)
    b = Velocity(2.0,2.0,3.0)
    @test a+b == Velocity(3.0,4.0,6.0)
    @test norm(b-a) == 1.0u"m/s"
    @test a*b == 15.0u"m^2/s^2"
    @test dot(a,b) == 15.0u"m^2/s^2"
    @test (a+1.0u"m/s")*(1.0u"m/s"+b) == 31.0u"m^2/s^2"
    @test (a-1.0u"m/s")*(1.0u"m/s"-b) == -5.0u"m^2/s^2"
    @test gety(a*2.0) == 4.0u"m/s"
    @test getx(a/2.0) == 0.5u"m/s"
    @test getz(2.0*a) == 6.0u"m/s"
    @test zero(a).y == 0.0u"m/s"

    @test normalize(Velocity(3.0,4.0,12.0)).x == 3.0/13.0*u"m/s"
    @test rotate_x(Velocity(0.0,1.0,0.0), 0.5pi).z == 1.0u"m/s"
    @test rotate_y(Velocity(0.0,0.0,1.0), 0.5pi).x == 1.0u"m/s"
    @test rotate_z(Velocity(1.0,0.0,0.0), 0.5pi).y == 1.0u"m/s"
    @test norm(pconvert([3.0,4.0,12.0],u"m/s")) == 13.0u"m/s"

    c = pconvert(array3d,u"m/s")
    @test mean(c) == Velocity(3.0,8.0,13.0)
    @test center(c) == Velocity(3.0,8.0,13.0)
    @test min_x(c) == 1.0u"m/s"
    @test max_x(c) == 5.0u"m/s"
    @test min_y(c) == 6.0u"m/s"
    @test max_y(c) == 10.0u"m/s"
    @test min_z(c) == 11.0u"m/s"
    @test max_z(c) == 15.0u"m/s"
end

@testset "3D Velocity Astro" begin
    a = VelocityAstro(1.0,2.0,3.0)
    b = VelocityAstro(2.0,2.0,3.0)
    @test a+b == VelocityAstro(3.0,4.0,6.0)
    @test norm(b-a) == 1.0u"kpc/Gyr"
    @test a*b == 15.0u"kpc^2/Gyr^2"
    @test dot(a,b) == 15.0u"kpc^2/Gyr^2"
    @test (a+1.0u"kpc/Gyr")*(1.0u"kpc/Gyr"+b) == 31.0u"kpc^2/Gyr^2"
    @test (a-1.0u"kpc/Gyr")*(1.0u"kpc/Gyr"-b) == -5.0u"kpc^2/Gyr^2"
    @test gety(a*2.0) == 4.0u"kpc/Gyr"
    @test getx(a/2.0) == 0.5u"kpc/Gyr"
    @test getz(2.0*a) == 6.0u"kpc/Gyr"
    @test zero(a).y == 0.0u"kpc/Gyr"

    @test normalize(VelocityAstro(3.0,4.0,12.0)).x == 3.0/13.0*u"kpc/Gyr"
    @test rotate_x(VelocityAstro(0.0,1.0,0.0), 0.5pi).z == 1.0u"kpc/Gyr"
    @test rotate_y(VelocityAstro(0.0,0.0,1.0), 0.5pi).x == 1.0u"kpc/Gyr"
    @test rotate_z(VelocityAstro(1.0,0.0,0.0), 0.5pi).y == 1.0u"kpc/Gyr"
    @test norm(pconvert([3.0,4.0,12.0],u"kpc/Gyr")) == 13.0u"kpc/Gyr"

    c = pconvert(array3d,u"kpc/Gyr")
    @test mean(c) == VelocityAstro(3.0,8.0,13.0)
    @test center(c) == VelocityAstro(3.0,8.0,13.0)
    @test min_x(c) == 1.0u"kpc/Gyr"
    @test max_x(c) == 5.0u"kpc/Gyr"
    @test min_y(c) == 6.0u"kpc/Gyr"
    @test max_y(c) == 10.0u"kpc/Gyr"
    @test min_z(c) == 11.0u"kpc/Gyr"
    @test max_z(c) == 15.0u"kpc/Gyr"
end

@testset "3D Acceleration" begin
    a = Acceleration(1.0,2.0,3.0)
    b = Acceleration(2.0,2.0,3.0)
    @test a+b == Acceleration(3.0,4.0,6.0)
    @test norm(b-a) == 1.0u"m/s^2"
    @test a*b == 15.0u"m^2/s^4"
    @test dot(a,b) == 15.0u"m^2/s^4"
    @test (a+1.0u"m/s^2")*(1.0u"m/s^2"+b) == 31.0u"m^2/s^4"
    @test (a-1.0u"m/s^2")*(1.0u"m/s^2"-b) == -5.0u"m^2/s^4"
    @test gety(a*2.0) == 4.0u"m/s^2"
    @test getx(a/2.0) == 0.5u"m/s^2"
    @test getz(2.0*a) == 6.0u"m/s^2"
    @test zero(a).y == 0.0u"m/s^2"

    @test normalize(Acceleration(3.0,4.0,12.0)).x == 3.0/13.0*u"m/s^2"
    @test rotate_x(Acceleration(0.0,1.0,0.0), 0.5pi).z == 1.0u"m/s^2"
    @test rotate_y(Acceleration(0.0,0.0,1.0), 0.5pi).x == 1.0u"m/s^2"
    @test rotate_z(Acceleration(1.0,0.0,0.0), 0.5pi).y == 1.0u"m/s^2"
    @test norm(pconvert([3.0,4.0,12.0],u"m/s^2")) == 13.0u"m/s^2"

    c = pconvert(array3d,u"m/s^2")
    @test mean(c) == Acceleration(3.0,8.0,13.0)
    @test center(c) == Acceleration(3.0,8.0,13.0)
    @test min_x(c) == 1.0u"m/s^2"
    @test max_x(c) == 5.0u"m/s^2"
    @test min_y(c) == 6.0u"m/s^2"
    @test max_y(c) == 10.0u"m/s^2"
    @test min_z(c) == 11.0u"m/s^2"
    @test max_z(c) == 15.0u"m/s^2"
end

@testset "3D Acceleration Astro" begin
    a = AccelerationAstro(1.0,2.0,3.0)
    b = AccelerationAstro(2.0,2.0,3.0)
    @test a+b == AccelerationAstro(3.0,4.0,6.0)
    @test norm(b-a) == 1.0u"kpc/Gyr^2"
    @test a*b == 15.0u"kpc^2/Gyr^4"
    @test dot(a,b) == 15.0u"kpc^2/Gyr^4"
    @test (a+1.0u"kpc/Gyr^2")*(1.0u"kpc/Gyr^2"+b) == 31.0u"kpc^2/Gyr^4"
    @test (a-1.0u"kpc/Gyr^2")*(1.0u"kpc/Gyr^2"-b) == -5.0u"kpc^2/Gyr^4"
    @test gety(a*2.0) == 4.0u"kpc/Gyr^2"
    @test getx(a/2.0) == 0.5u"kpc/Gyr^2"
    @test getz(2.0*a) == 6.0u"kpc/Gyr^2"
    @test zero(a).y == 0.0u"kpc/Gyr^2"

    @test normalize(AccelerationAstro(3.0,4.0,12.0)).x == 3.0/13.0*u"kpc/Gyr^2"
    @test rotate_x(AccelerationAstro(0.0,1.0,0.0), 0.5pi).z == 1.0u"kpc/Gyr^2"
    @test rotate_y(AccelerationAstro(0.0,0.0,1.0), 0.5pi).x == 1.0u"kpc/Gyr^2"
    @test rotate_z(AccelerationAstro(1.0,0.0,0.0), 0.5pi).y == 1.0u"kpc/Gyr^2"
    @test norm(pconvert([3.0,4.0,12.0],u"kpc/Gyr^2")) == 13.0u"kpc/Gyr^2"

    c = pconvert(array3d,u"kpc/Gyr^2")
    @test mean(c) == AccelerationAstro(3.0,8.0,13.0)
    @test center(c) == AccelerationAstro(3.0,8.0,13.0)
    @test min_x(c) == 1.0u"kpc/Gyr^2"
    @test max_x(c) == 5.0u"kpc/Gyr^2"
    @test min_y(c) == 6.0u"kpc/Gyr^2"
    @test max_y(c) == 10.0u"kpc/Gyr^2"
    @test min_z(c) == 11.0u"kpc/Gyr^2"
    @test max_z(c) == 15.0u"kpc/Gyr^2"
end
