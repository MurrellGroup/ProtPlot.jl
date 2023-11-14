@testset "segment.jl" begin

    using .Ribbon

    @testset "Segment" begin
        chain = Chain("A", Backbone(randn(3, 4, 5)))
        chain.ssvector .= [Loop, Loop, Loop, Helix, Strand]
        segment = Segment{Loop}(chain, 1:3)
        @test segment.chain == chain
        @test segment.range == 1:3
    end

    @testset "extend_segment" begin
        chain = Chain("A", Backbone(randn(3, 4, 5)))
        chain.ssvector .= [Loop, Helix, Strand, Strand, Helix]
        segment = Segment{Strand}(chain, 3:4)
        @test extend_segment(segment, 0:3) == Segment{MiSSing}(chain, 2:5)
    end

    @testset "segments" begin
        coords = randn(3, 4, 3)
        chain = Chain("B", Backbone(coords))
        chain.ssvector .= [Loop, Helix, Helix]
        chain_segments = segments(chain)
        @test length(chain_segments) == 2
        @test chain_segments[1] isa Segment{Loop, Float64}
        @test chain_segments[1].chain == chain
        @test chain_segments[1].range == 1:1
        @test chain_segments[1].backbone.coords == coords[:, :, 1:1]
        @test chain_segments[2].backbone.coords == coords[:, :, 2:3]
    end

end