import GLMakie

@testset "render.jl" begin
    protein = readpdb("data/1ZAK.pdb")
    @test ribbon(protein) isa GLMakie.Scene
end