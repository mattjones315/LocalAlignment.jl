function build_arg_table()
    settings = ArgParseSettings(
        description="\033[32mLocal Sequence Alignment\033[0m\n\n\n\n" *
        "\033[31m Dynamic Programming algorithms for local sequence alignment\033[0m"
    )

    @add_arg_table settings begin
        "align"
            action = :command
            help = "run package from command line arguments"
        "ls"
            action = :command
            help = "list all modified simulations available"
    end

    settings["ls"].description = "Prints out the chosen algorithm for sequence alginment. "

    @add_arg_table settings["align"] begin
        "--addprocs", "-p"
            help = "Add additional processors"
            arg_type = Int
            default = 0
        "--algorithm", "-a"
            help = "Which algorithm to run"
            arg_type = Compat.String
            default = "smithwaterman"
        "--gap_penalty", "-g"
            help = "Gap penalty"
            arg_type = Float64
            default = 1.0
        "--extension_penalty", "-e"
            help = "Penalty for extending a gap"
            arg_type = Float64
            default = 1.0
        "--score_mat", "-s"
            help = "File path for scoring matrix"
            required = true
        "sequence_1"
            help = "Fasta file 1 for sequence alignemnt"
            required = true
        "sequence_2"
            help = "Fasta file 2 for sequence alignment"
            required = true
        "output_file"
            help = "File to output results to [.CSV, .TSV, etc]"
            required = true
    end

    settings["align"].description

    if typeof(Base.source_dir()) != Void
        settings.epilog = readstring(normpath(joinpath(Base.source_dir(),"..","LICENSE")))
    end

    return settings
end
