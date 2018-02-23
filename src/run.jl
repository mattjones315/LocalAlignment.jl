using ArgParse
using Compat
import Compat: UTF8String, ASCIIString, view, readstring, foreach

include("CLI.jl")
include("commands.jl")

function main()
    parsed_args = parse_args(build_arg_table())
    command = parsed_args["%COMMAND%"]

    if command == "help"
        foreach(x -> println(x), help())
    elseif command == "align_many"
        align_many(parsed_args[command]["algorithm"],
        parsed_args[command]["gap_penalty"],
        parsed_args[command]["extension_penalty"],
        parsed_args[command]["sequences"],
        parsed_args[command]["score_mat"],
        parsed_args[command]["output_file"],
        parsed_args[command]["write_alignments"]
        )
    elseif command == "compare_pr"
        compare_pr(parsed_args[command]["algorithm"],
        parsed_args[command]["pos_seqs"],
        parsed_args[command]["neg_seqs"],
        parsed_args[command]["score_mat"],
        parsed_args[command]["gp_min"],
        parsed_args[command]["gp_max"],
        parsed_args[command]["output_file"],
        parsed_args[command]["tpr"])
    elseif command == "roc"
        bench_roc(parsed_args[command]["algorithm"],
        parsed_args[command]["pos_seqs"],
        parsed_args[command]["neg_seqs"],
        parsed_args[command]["score_mat"],
        parsed_args[command]["gap_penalty"],
        parsed_args[command]["extension_penalty"],
        parsed_args[command]["output_file"])
    else
        align(parsed_args[command]["algorithm"],
        parsed_args[command]["gap_penalty"],
        parsed_args[command]["extension_penalty"],
        parsed_args[command]["sequence_1"],
        parsed_args[command]["sequence_2"],
        parsed_args[command]["score_mat"],
        parsed_args[command]["output_file"]
        )
    end
end

# fire up simulation if run using command line
if !isinteractive()
    main()
end
