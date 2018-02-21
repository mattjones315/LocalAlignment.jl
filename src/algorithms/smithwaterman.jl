using DataFrames;

"""
sw_align(seq1::String, seq2::String, smat::DataFrame, gapp::Float64, gapcp::Float64, out_fp::String)
Align FASTA sequences SEQ1 and SEQ2 using the Smith-Waterman algorithm with gap penalty
GP and extension penalty EP and according to the scoring matrix SMAT.
Stores the result in OUT_FP.

"""
function sw_align(seq1, seq2, smat, gp, ep, out_fp; write_align = true)

    M = length(seq1);
    N = length(seq2);

    H = zeros(M+1, N+1);
    Wk = zeros(M+1);
    Wm = zeros(N+1);

    for i in 2:(M+1)
        for j in 2:(N+1)

            diag = (H[i-1, j-1] + smat[smat[:, :AA] .== string(uppercase(seq1[i-1])),
                                            Symbol(uppercase(seq2[j-1]))])[1];

            # W_k = extension penalty * length of extension + gap penalty
            # We don't take into account extension length until after a gap is allowed
            vert = maximum([(H[i-k, j] - (k)*ep - gp) for k in 1:(i-1)]);
            horiz = maximum([(H[i, j-k] - (k)*ep - gp) for k in 1:(j-1)]);
            opt = max(diag, vert, horiz, 0);
            H[i, j] = opt;
        end
    end

    #println(DataFrame(H));
    if write_align
        sw_traceback(seq1, seq2, H, out_fp)
    else
        return findmax(H)[1]
    end

end

function many_sw_align(seqs, smat, gp, ep, out_fp)

    scores = [];

    progress_i = 1

    for seq_pair in seqs

        println(string("Aligning pair ", progress_i, " of ", length(seqs)))

        s = sw_align(seq_pair[1], seq_pair[2], smat, gp, ep, out_fp; write_align=false)
        push!(scores, s)

        progress_i += 1

    end

    open(out_fp, "w") do f
        for s in scores
            write(f, string(s, '\t'));
        end
        write(f, '\n');
    end

end

function sw_traceback(seq1, seq2, H, out_fp)

    m = length(seq1)+1;
    n = length(seq2)+1;

    align1 = "";
    matches = "";
    align2 = "";

    mii = ind2sub(H, findmax(H)[2]);

    i, j = mii[1], mii[2];

    while m > i
        align2 = string("-", align2)
        matches = string(" ", matches);
        align1 = string(seq1[m-1], align1);
        m -= 1
    end

    while n > j
        align1 = string("-", align1)
        matches = string(" ", matches);
        align2 = string(seq2[n-1], align2);
        n -= 1
    end

    while H[i, j] != 0

        diag = H[i-1, j-1];
        up = H[i-1, j];
        left = H[i, j-1];

        if diag >= up && diag >= left
            matches = string("|", matches);
            align1 = string(seq1[i-1], align1);
            align2 = string(seq2[j-1], align2);
            i -= 1;
            j -= 1;
        elseif up >= left
            matches = string(" ", matches);
            align1 = string(seq1[i-1], align1);
            align2 = string("-", align2)
            i -= 1
        else
            matches = string(matches, " ");
            align1 = string(align1, seq1[i-1]);
            align2 = string(align2, "-")
            j -= 1
        end
    end

    open(out_fp, "w") do f
        #write(f, string(seq1, " aligned to ", seq2, '\n'))
        write(f, string(align1, '\n'))
        write(f, string(matches, '\n'))
        write(f, string(align2, '\n'))
        write(f, string("Score: ", findmax(H)[1], '\n'))
    end

    return

end
