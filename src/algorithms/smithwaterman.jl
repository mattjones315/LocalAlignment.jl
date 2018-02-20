using DataFrames;

"""
sw_align(seq1::String, seq2::String, smat::DataFrame, gapp::Float64, gapcp::Float64, out_fp::String)
Align FASTA sequences SEQ1 and SEQ2 using the Smith-Waterman algorithm with gap penalty
GP and extension penalty EP and according to the scoring matrix SMAT.
Stores the result in OUT_FP.

"""
function sw_align(seq1, seq2, smat, gp, ep, out_fp)

    M = length(seq1);
    N = length(seq2);

    H = zeros(M+1, N+1);
    for i in 2:(M+1)
        for j in 2:(N+1)

            diag = (H[i-1, j-1] + smat[smat[:, :AA] .== string(seq1[i-1]), Symbol(seq2[j-1])])[1];

            # W_k = extension penalty * length of extension + gap penalty
            # We don't take into account extension length until after a gap is allowed
            vert = max.([(H[i-k, j] - (k-1)*ep - gp) for k in 1:(i-1)])[1];
            horiz = max.([(H[i, j-m] - (m-1)*ep - gp) for m in 1:(j-1)])[1];
            opt = max.(diag, vert, horiz, 0);
            H[i, j] = opt[1];
        end
    end

    sw_traceback(seq1, seq2, H, out_fp)

end

function sw_traceback(seq1, seq2, H, out_fp)

    N = length(seq1);
    M = length(seq2);

    i, j = N, M;

    while i > 0 or j > 0
        println(i);
    end


    return

end
