@testset "segment.jl" begin

    using .Ribbon

    @testset "Segment" begin
        chain = Protein.Chain("A", Backbone(randn(3, 15)))
        chain.ssvector .= ['-', '-', '-', 'H', 'E']
        segment = Segment{Ribbon.Coil}(chain, 1:3)
        @test segment.chain == chain
        @test segment.range == 1:3
    end

    @testset "extend_segment" begin
        chain = Protein.Chain("A", Backbone(randn(3, 15)))
        chain.ssvector .= ['-', 'H', 'E', 'E', 'H']
        segment = Segment{Ribbon.Strand}(chain, 3:4)
        @test extend_segment(segment, 0:3) == Segment{Ribbon.Unassigned}(chain, 2:5)
    end

    @testset "segments" begin
        coords = randn(3, 9)
        chain = Protein.Chain("B", Backbone(coords))
        chain.ssvector .= ['-', 'H', 'H']
        chain_segments = segments(chain)
        @test length(chain_segments) == 2
        @test chain_segments[1] isa Segment{Ribbon.Coil}
        @test chain_segments[1].chain == chain
        @test chain_segments[1].range == 1:1
        @test chain_segments[1].backbone.coords == coords[:, 1:3]
        @test chain_segments[2].backbone.coords == coords[:, 4:9]
    end

end