import GLMakie

@testset "render.jl" begin
    protein = readpdb("data/1ZAK.pdb")
    @test ribbon(protein) isa GLMakie.Scene
    scene = ribbon(protein; camcontrols=(; lookat=Vec3f(30, 0, 60), eyeposition=Vec3f(160, -75, 0), upvector=Vec3f(0, 0, 1)))
    @test scene.camera.eyeposition[] == Vec3f(160, -75, 0)
    @test scene.camera.lookat[] == Vec3f(30, 0, 60)
end
