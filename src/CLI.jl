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
        "align_many"
            action = :command
            help = "run package for multiple different sequences at once"
        "compare_pr"
            action = :command
            help = "run package for benchmarking the false positive rate"
        "roc"
            action = :command
            help = "run package for generating ROC data between pos and neg pairs"
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

    @add_arg_table settings["align_many"] begin
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
        "--write_alignments"
            help = "Write all alignments"
            default = false
        "sequences"
            help = "File containing list of FASTA files"
            required = true
        "output_file"
            help = "File to output results to [.CSV, .TSV, etc]"
            required = true
    end

    settings["align_many"].description

    @add_arg_table settings["compare_pr"] begin
        "--addprocs", "-p"
            help = "Add additional processors"
            arg_type = Int
            default = 0
        "--algorithm", "-a"
            help = "Which algorithm to run"
            arg_type = Compat.String
            default = "smithwaterman"
        "--score_mat", "-s"
            help = "File path for scoring matrix"
            required = true
        "--tpr"
            help = "true positive rate for thresholding"
            arg_type = Float64
            default = 0.7
        "--gp_min"
            help = "minimum gap penalty to test"
            arg_type = Int64
            default = 1
        "--gp_max"
            help = "maximum gap penalty to test"
            arg_type = Int64
            default = 20
        "pos_seqs"
            help = "File containing list of 'positive' FASTA files"
            required = true
        "neg_seqs"
            help = "File containing a list of 'negative' FASTA files"
            required = true
        "output_file"
            help = "File to output results to [.CSV, .TSV, etc]"
            required = true
    end

    settings["compare_pr"].description

    @add_arg_table settings["roc"] begin
        "--addprocs", "-p"
            help = "Add additional processors"
            arg_type = Int
            default = 0
        "--algorithm", "-a"
            help = "Which algorithm to run"
            arg_type = Compat.String
            default = "smithwaterman"
        "--score_mat", "-s"
            help = "File path for scoring matrix"
            required = true
        "--gap_penalty", "-g"
            help = "Gap penalty"
            arg_type = Float64
            default = 1.0
        "--extension_penalty", "-e"
            help = "Penalty for extending a gap"
            arg_type = Float64
            default = 1.0
        "pos_seqs"
            help = "File containing list of 'positive' FASTA files"
            required = true
        "neg_seqs"
            help = "File containing a list of 'negative' FASTA files"
            required = true
        "output_file"
            help = "File to output results to [.CSV, .TSV, etc]"
            required = true
    end

    settings["roc"].description

    if typeof(Base.source_dir()) != Void
        settings.epilog = readstring(normpath(joinpath(Base.source_dir(),"..","LICENSE")))
    end

    return settings
end
