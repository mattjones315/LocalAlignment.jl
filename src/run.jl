using ArgParse
using YAML
using Compat
import Compat: UTF8String, ASCIIString, view, readstring, foreach

include("CLI.jl")
include("commands.jl")

function main()
    parsed_args = parse_args(build_arg_table())
    command = parsed_args["%COMMAND%"]

    if command == "ls"
        foreach(x -> println(x), ls())
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