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
            vert = maximum([(H[i-k, j] - (k-1)*ep - gp) for k in 1:(i-1)]);
            horiz = maximum([(H[i, j-k] - (k-1)*ep - gp) for k in 1:(j-1)]);
            opt = max(diag, vert, horiz, 0);
            H[i, j] = opt;
        end
    end

    #println(DataFrame(H));
    if write_align
        sw_traceback(seq1, seq2, H, out_fp)
    else
        return findmax(H)[1], H
    end

end

function many_sw_align(seqs, smat, gp, ep, out_fp; write_output=true, normalize=false,
                        write_alignments = true)

    println(gp)
    println(ep)

    scores = [];
    Hs = [];

    progress_i = 1

    for seq_pair in seqs

        println(string("Aligning pair ", progress_i, " of ", length(seqs)))

        s, H = sw_align(seq_pair[1], seq_pair[2], smat, gp, ep, out_fp; write_align=false)

        push!(Hs, H)

        if normalize
            s = s / min(length(seq_pair[1]), length(seq_pair[2]))
        end

        push!(scores, s)

        progress_i += 1

    end

    if write_output && !write_alignments
        open(out_fp, "w") do f
            for s in scores
                write(f, string(s, '\t'));
            end
            write(f, '\n');
        end
    elseif write_alignments && write_output
        sw_traceback_many(seqs, Hs, out_fp)
    else
        scores
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

    _end = false
    while !_end
        diag = H[i-1, j-1];
        up = H[i-1, j];
        left = H[i, j-1];

        if diag >= up && diag >= left
            if diag == 0
                _end = true
            else
                align1 = string(seq1[i-1], align1)
                align2 = string(seq2[j-1], align2)
                i -= 1
                j -= 1
            end
        elseif up > diag && up >= left
            if up == 0
                _end = true
            else
                align1 = string(seq1[i-1], align1)
                align2 = string("-", align2)
                i -= 1
            end
        else
            if left == 0
                _end = true
            else
                align1 = string("-", align1)
                align2 = string(seq2[j-1], align2)
                j -= 1
            end
        end
    end

    align1 = string(seq1[i-1], align1)
    align2 = string(seq2[j-1], align2)

    open(out_fp, "w") do f
        write(f, string(align1, '\n'))
        write(f, string(align2, '\n'))
        write(f, string("Score: ", findmax(H)[1], '\n'))
    end

    return

end

function sw_traceback_many(seqs, Hs, out_fp)

    open(out_fp, "w") do f
        for i in 1:length(seqs)

            seq1 = seqs[i][1];
            seq2 = seqs[i][2];
            H = Hs[i];

            m = length(seq1)+1;
            n = length(seq2)+1;

            align1 = "";
            matches = "";
            align2 = "";

            mii = ind2sub(H, findmax(H)[2]);

            i, j = mii[1], mii[2];

            _end = false
            while !_end
                diag = H[i-1, j-1];
                up = H[i-1, j];
                left = H[i, j-1];

                if diag >= up && diag >= left
                    if diag == 0
                        _end = true
                    else
                        align1 = string(seq1[i-1], align1)
                        align2 = string(seq2[j-1], align2)
                        i -= 1
                        j -= 1
                    end
                elseif up > diag && up >= left
                    if up == 0
                        _end = true
                    else
                        align1 = string(seq1[i-1], align1)
                        align2 = string("-", align2)
                        i -= 1
                    end
                else
                    if left == 0
                        _end = true
                    else
                        align1 = string("-", align1)
                        align2 = string(seq2[j-1], align2)
                        j -= 1
                    end
                end
            end

            align1 = string(seq1[i-1], align1)
            align2 = string(seq2[j-1], align2)

            write(f, string(align1, '\n'))
            write(f, string(align2, '\n'))
            write(f, string("Score: ", findmax(H)[1], '\n'))

        end
    end

end
