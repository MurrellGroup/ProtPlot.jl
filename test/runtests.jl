using GLMakie
using ProtPlot
using Test

using ProteinChains
import BioStructures

@testset "ProtPlot.jl" begin
    @test ribbon(pdb"1ZAK") isa Makie.FigureAxisPlot

    # gap in chain
    @test ribbon(pdb"1M4X") isa Makie.FigureAxisPlot

    # single-residue subchain
    @test ribbon(pdb"1M4X"A[[1:10; 20; 30:40]]) isa Makie.FigureAxisPlot

    camcontrols = (; lookat=Vec3f(30, 0, 60), eyeposition=Vec3f(160, -75, 0), upvector=Vec3f(0, 0, 1))
    scene = ribbon_scene(pdb"1ASS"; camcontrols)
    @test scene.camera.eyeposition[] == Vec3f(160, -75, 0)

    @test atomplot(pdb"1ASS") isa Makie.FigureAxisPlot
    @test atomplot(pdb"1ASS"A) isa Makie.FigureAxisPlot
    @test atomplot(Frames(pdb"1ASS"A)) isa Makie.FigureAxisPlot

    @test ramachandran(pdb"1ASS") isa Makie.FigureAxisPlot

    mktempdir() do dir
        BioStructures.downloadpdb("1M4X"; dir)
        @test atomplot(read(joinpath(dir, "1M4X.pdb"), BioStructures.PDBFormat)) isa Makie.FigureAxisPlot
    end
end
