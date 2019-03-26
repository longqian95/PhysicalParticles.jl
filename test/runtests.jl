include("../src/PhysicalParticles.jl")

using .PhysicalParticles
using Test, Unitful, UnitfulAstro

array2d = [1.0 2.0 3.0 4.0 5.0;
           6.0 7.0 8.0 9.0 10.0]
array3d = [1.0 2.0 3.0 4.0 5.0;
           6.0 7.0 8.0 9.0 10.0;
           11.0 12.0 13.0 14.0 15.0]

@testset "Normal Points" begin
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
end

@testset "2D Physical Vectors" begin

end

@testset "3D Physical Vectors" begin

end

@testset "Physical Particles" begin

end
