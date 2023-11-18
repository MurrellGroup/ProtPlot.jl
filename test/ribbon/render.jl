import GLMakie

@testset "render.jl" begin
    protein = pdb_to_protein("data/1ZAK.pdb")
    @test ribbon(protein) isa GLMakie.Scene
end