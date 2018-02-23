using DataFrames;
using Distributions;

#source_dir = typeof(Base.source_dir()) == Void ? joinpath(Pkg.dir("LocalAlignment"), "src") : Base.source_dir()
source_dir = pwd()

include(joinpath(source_dir, "src", "algorithms", "smithwaterman.jl"));
include(joinpath(source_dir, "src", "utils", "parse.jl"));
include(joinpath(source_dir, "src", "utils", "benchmarking.jl"));

function optimize(S, T, N, g, e, pos_a_fp, neg_a_fp)

    pos_aligns = parse_alignments(pos_a_fp);
    neg_aligns = parse_alignments(neg_a_fp);

    Scurr = S

    prev_pos_scores = []
    prev_neg_scores = []
    for pos in pos_aligns
        push!(prev_pos_scores, score_alignment(pos[1], pos[2], Scurr, g, e)[1])
    end

    for neg in neg_aligns
        push!(prev_neg_scores, score_alignment(neg[1], neg[2], Scurr, g, e)[1])
    end

    Jprev = eval_objective_function(prev_pos_scores, prev_neg_scores)


    for i in 1:N
        println(string("Iteration: ", i));
        curr_pos_scores = []
        curr_neg_scores = []

        Sp = gen_new_matrix(Scurr)

        for pos in pos_aligns
            push!(curr_pos_scores, score_alignment(pos[1], pos[2], Sp, g, e)[1])
        end

        for neg in neg_aligns
            push!(curr_neg_scores, score_alignment(neg[1], neg[2], Sp, g, e)[1])
        end

        Jcurr = eval_objective_function(curr_pos_scores, curr_neg_scores)

        if Jcurr >= Jprev
            Scurr = Sp
            prev_neg_scores = curr_neg_scores
            prev_pos_scores = curr_pos_scores
            Jprev = Jcurr
        else
            p = e^( (Jcurr - Jprev) * i ^T)
            if rand(Uniform()) < p
                Scurr = Sp
                prev_neg_scores = curr_neg_scores
                prev_pos_scores = curr_pos_scores
                Jprev = Jcurr
            else
                continue;
            end
        end

    end

    Scurr

end

function eval_objective_function(pos_scores, neg_scores)

    val = 0.0
    fpr = [0.0, 0.1, 0.2, 0.3]

    for fp in fpr
        thresh = find_tp_threshold(neg_scores, fp)
        val += length(pos_scores[pos_scores .> thresh]) / length(pos_scores)
    end

    return val
end

function score_alignment(align1, align2, S, g, e)

    score = 0.0
    penalty = 0.0

    for i in 1:length(align1)

        if align1[i] != '-' && align2[i] != '-'
            score += S[S[:, :AA] .== string(uppercase(align1[i])),
                            Symbol(uppercase(align2[i]))];
        end

    end

    curr_e = 0
    for i in 1:length(align1)

        if align1[i] == '-' || align2[i] == '-'
            if curr_e == 0
                penalty += g
            else
                penalty += e
            end
            curr_e += 1
        else
            curr_e = 0;
        end

    end

    score = score - penalty
    return score

end

function parse_alignments(fp)

    alignments = []

    f = open(fp, "r");
    lines = collect(readlines(f))

    for i in 1:3:length(lines)

        align1 = strip(lines[i]);
        align2 = strip(lines[i+1]);
        push!(alignments, [align1, align2]);

    end

    alignments

end

function gen_new_matrix(S)

    Sp = deepcopy(S)
    Sp = Array(Sp[:,2:end])
    aa = S[:,:AA]

    for i in 1:size(Sp, 1)
        for j in size(Sp, 1):-1:i
            Sp[i, j] += randn()
            Sp[j, i] = Sp[i,j]
        end
    end

    Sp = DataFrame(Sp);
    Sp = hcat(aa, Sp);
    names!(Sp, names(S));

    Sp

end
