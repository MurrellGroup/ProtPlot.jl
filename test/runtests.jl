using Makie
using ProtPlot
using Test

using Backboner

@testset "ProtPlot.jl" begin
    chains = readpdb("data/1ZAK.pdb")
    @test ribbon(chains) isa Makie.AbstractPlot # bad test but at least it runs

    scene = ribbon(protein; camcontrols=(; lookat=Vec3f(30, 0, 60), eyeposition=Vec3f(160, -75, 0), upvector=Vec3f(0, 0, 1)))
    @test scene.camera.eyeposition[] == Vec3f(160, -75, 0)
end
