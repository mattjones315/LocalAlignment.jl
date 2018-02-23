load_file1 = joinpath("src", "algorithms", "smithwaterman.jl")
include(normpath(joinpath(Base.source_dir(), "..", "src", "utils", "benchmarking")))
include(normpath(joinpath(Base.source_dir(),"..",load_file)))

using Base.Test
using DataStructures
using ColorBrewer
using Gadfly
using Compat
import Compat: UTF8String, ASCIIString, view

println("Running tests:")

tic()
println("Running Smith-Waterman Tests")
@test include("test_sm.jl")
println("Running Benchmarking Tests")
@test include("test_benchmarking.jl")
toc()
