using GLMakie
using ProtPlot
using Test

@testset "ProtPlot.jl" begin
    @test ribbon_scene(pdb"1ZAK") isa Scene # shitty test but at least it confirms it doesn't throw an error

    camcontrols = (; lookat=Vec3f(30, 0, 60), eyeposition=Vec3f(160, -75, 0), upvector=Vec3f(0, 0, 1))
    scene = ribbon_scene(pdb"1ASS"; camcontrols)
    @test scene.camera.eyeposition[] == Vec3f(160, -75, 0)
end
