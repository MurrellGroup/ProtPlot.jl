using GLMakie
using ProtPlot
using Test

@testset "ProtPlot.jl" begin
    @test ribbon(pdb"1ZAK") isa Makie.FigureAxisPlot

    camcontrols = (; lookat=Vec3f(30, 0, 60), eyeposition=Vec3f(160, -75, 0), upvector=Vec3f(0, 0, 1))
    scene = ribbon_scene(pdb"1ASS"; camcontrols)
    @test scene.camera.eyeposition[] == Vec3f(160, -75, 0)

    @test atomplot(pdb"1ASS") isa Makie.FigureAxisPlot

    @test ramachandran(pdb"1ASS") isa Makie.FigureAxisPlot
end
