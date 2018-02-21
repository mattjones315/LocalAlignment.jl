using DataFrames;
using CSV;

function benchmark_fpr_vs_tpr(pos_seqs, neg_seqs, smat, out_fp, tpr)

    gps = [i for i in 1:20];
    eps = [i for i in 1:5];

    G = length(gps);
    E = length(eps);

    fpr_df = zeros(G, E);

    for g in 1:G
        for e in 1:E

            gp = gps[g];
            ep = eps[e];
            println(string("Gap Penalty: ", gp, " Extension Penalty: ", ep))
            pos_scores = many_sw_align(pos_seqs, smat, gp, ep, out_fp; write_output=false);
            neg_scores = many_sw_align(neg_seqs, smat, gp, ep, out_fp; write_output=false);

            thresh = find_tp_threshold(pos_scores, tpr);
            println(thresh)
            fpr = length(neg_scores[neg_scores .> thresh]) / length(neg_scores)

            fpr_df[g, e] = fpr

        end
    end
    fpr_df = DataFrame(fpr_df)
    names!(fpr_df, [Symbol("$i") for i in eps])
    gpdf = DataFrame(GP = gps)
    fpr_df = hcat(gpdf, fpr_df);
    CSV.write(out_fp, fpr_df, delim='\t');
end

function find_tp_threshold(scores, tpr)

    for i in minimum(scores):maximum(scores)

        rate = length(scores[scores .> i]) / length(scores)
        if rate <= tpr
            return i
        end
    end

    return i

end
