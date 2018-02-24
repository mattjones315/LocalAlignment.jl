using DataFrames;
using CSV;

"""
benchmark_fpr_vs_trp(pos_seqs::Array, neg_seqs::Array, smat::DataFrame, gp_min::Float64,
                        gp_max::Float64, out_fp::String, tpr::Float64)
Benchmark the False Positive Rates at the given true positive rate TPR for positive
sequences in POS_SEQS and negative sequences in NEG_SEQS subject to substitution
matrix SMAT. Iterate over parameters -- extension penalties between 1 and 5, and
gap penalties between GP_MIN and GP_MAX. Write a table out to OUT_FP storing
the false positive rates for a given set of parameters.
"""

function benchmark_fpr_vs_tpr(pos_seqs, neg_seqs, smat, gp_min, gp_max, out_fp, tpr)

    gps = [i for i in gp_min:gp_max];
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

"""
find_tp_threshold(scores::Array, tpr::Float64)
Find the score threshold in SCORES for a true positive rate of TPR.
"""
function find_tp_threshold(scores, tpr)
    for i in linspace(minimum(scores), maximum(scores), 200)

        rate = length(scores[scores .> i]) / length(scores)
        if rate <= tpr
            return i
        end
    end

    return i

end

"""
roc(pos_seqs::Array, neg_seqs::Array, smat::DataFrame, gp::Float64, ep::Float64, out_fp::String)
For a given gap penalty GP, extension penalty EP, and scoring matrix SMAT, construct the
ROC curve for pos_seqs and neg_seqs. Write output to OUT_FP consisting of coordinates in the form
(FPR, TPR).
"""
function roc(pos_seqs, neg_seqs, smat, gp, ep, out_fp)

    pos_scores = many_sw_align(pos_seqs, smat, gp, ep, out_fp; write_output=false)
    neg_scores = many_sw_align(neg_seqs, smat, gp, ep, out_fp; write_output=false)

    tprs = collect(linspace(0, 1, 50));

    roc_df = zeros(length(tprs), 2);

    for i in 1:length(tprs)
        thresh = find_tp_threshold(pos_scores, tprs[i]);
        fpr = length(neg_scores[neg_scores .> thresh]) / length(neg_scores)
        roc_df[i, :] = [fpr, tprs[i]];

    end

    if maximum(roc_df[:,1]) < 1
        extra_space = collect(linspace(maximum(roc_df[:,1]), 1.0, 10))
        for e in 1:length(extra_space)
            roc_df = vcat(roc_df, [extra_space[e] 1.0]);
        end
    end

    open(out_fp, "w") do f
        write(f, "FPR\tTPR\n");
        for i in 1:size(roc_df, 1)
            write(f, string(roc_df[i,1], '\t', roc_df[i,2], '\n'));
        end
    end

end
