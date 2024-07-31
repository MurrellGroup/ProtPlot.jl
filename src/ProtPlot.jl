module ProtPlot

export Ribbon, ribbon, ribbon!
export ribbon_scene

export Ramachandran, ramachandran, ramachandran!

export PDBEntry

using Makie
using ColorTypes

import Backboner
import Backboner.Protein: readpdb
export Backboner, readpdb

import BioStructures

include("Ribbon/Ribbon.jl")
include("ramachandran.jl")

PDBEntry(pdbid::AbstractString; format=BioStructures.PDBFormat, kwargs...) = mktempdir() do dir
    path = BioStructures.downloadpdb(pdbid; dir, format, kwargs...)
    Backboner.Protein.readchains(path, format)
end

end
