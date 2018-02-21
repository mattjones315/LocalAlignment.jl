# Sane behavior when run from the REPL
source_dir = typeof(Base.source_dir()) == Void ? joinpath(Pkg.dir("LocalAlignment"), "src") : Base.source_dir()

"""
align(algorithm::Compat.String, gap_penalty::Float64, gap_cont_penalty::Float64,
            seq1::Compat.String, seq2::Compat.String, smat::Compat.String, outputfile::Compat.String)
Align SEQ1 to SEQ2 using the scoring matrix SMAT and paramters GAP_PENALTY and GAP_CONT_PENALTY
for the ALGORITHM defined by user.
"""
function align(algorithm::Compat.String,
                gap_penalty::Float64,
                extension_penalty::Float64,
                seq1::Compat.String,
                seq2::Compat.String,
                smat::Compat.String,
                outputfile::Compat.String)

    if algorithm == "smithwaterman"
        include(joinpath(source_dir, "algorithms", "smithwaterman.jl"));
        include(joinpath(source_dir, "utils", "parse.jl"));
        # run
        Base.invokelatest(sw_align, Base.invokelatest(parse_input, seq1, seq2, smat)...,
                        gap_penalty, extension_penalty, outputfile);
    end

end

"""
align(algorithm::Compat.String, gap_penalty::Float64, gap_cont_penalty::Float64,
            seq1::Compat.String, seq2::Compat.String, smat::Compat.String, outputfile::Compat.String)
Align SEQ1 to SEQ2 using the scoring matrix SMAT and paramters GAP_PENALTY and GAP_CONT_PENALTY
for the ALGORITHM defined by user.
"""
function align_many(algorithm::Compat.String,
                gap_penalty::Float64,
                extension_penalty::Float64,
                seqs::Compat.String,
                smat::Compat.String,
                outputfile::Compat.String)

    if algorithm == "smithwaterman"
        include(joinpath(source_dir, "algorithms", "smithwaterman.jl"));
        include(joinpath(source_dir, "utils", "parse.jl"));
        # run
        Base.invokelatest(many_sw_align, Base.invokelatest(parse_seq_list, seqs, smat)...,
                        gap_penalty, extension_penalty, outputfile);
    end

end
