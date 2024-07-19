module ProtPlot

using Makie
using ColorTypes

import Backboner
import Backboner.Protein: readpdb
export Backboner, readpdb

import BioStructures

function PDBEntry(pdbid::AbstractString; format=BioStructures.PDBFormat, kwargs...)
    chains = mktempdir() do temp_dir
        path = BioStructures.downloadpdb(pdbid; dir=temp_dir, kwargs...)
        Backboner.Protein.readchains(path, format)
    end
    return chains
end

export PDBEntry

include("Ribbon/Ribbon.jl")
include("ramachandran.jl")

end
