load_file1 = joinpath("src", "algorithms", "smithwaterman.jl")
source_dir = pwd()

include(joinpath(source_dir, "../src", "algorithms", "smithwaterman.jl"));
include(joinpath(source_dir, "../src", "utils", "parse.jl"));
include(joinpath(source_dir, "../src", "utils", "benchmarking.jl"));

# include(normpath(joinpath(Base.source_dir(), "..", "src", "utils", "benchmarking.jl")))
# include(normpath(joinpath(Base.source_dir(),"..",load_file)))

using Base.Test
using Compat
import Compat: UTF8String, ASCIIString, view

function test_tpr()

    scores = [1,2,3,4,5,6,7,8,9,10]
    thresh = find_tp_threshold(scores, .7)

    @test thresh < 4

    scores = [.1, .2, .3, .4, .5, .6, .7, 8., .9, .10]
    tresh = find_tp_threshold(scores, .7)

    @test thresh < .4

    return

end

function test_sw_score()

    seq1 = "PRTEINS";
    seq2 = "PRTWPSEIN";

    smat = parse_parse_score_matrix(joinpath(source_dir, "../scoring", "BLOSUM50"))

    score = sw_align(seq1, seq2, smat, 1, 1, "test.txt"; write_align = false);

    @test score = 37.0


    score = sw_align(seq1, seq2, smat, 5, 3, "test.txt"; write_align = false);

    @test score = 29.0

    return



end



println("Running tests:")
test_tpr();
