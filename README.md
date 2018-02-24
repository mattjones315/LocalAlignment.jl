# bmi203-hw-3

| **Build Status** |
|:---:|
| [![][travis-img]][travis-url] |

Dynamic Programming local alignment using the Smith-Waterman algorithm. Takes in FASTA sequence files and will return an optimal alignment of the the sequences.

[travis-img]: http://img.shields.io/travis/mattjones315/bmi203-hw-3.svg
[travis-url]: https://travis-ci.org/mattjones315/bmi203-hw-3

# Getting Started

You can clone this repository with

```
git clone https://github.com/mattjones315/bmi203-hw-3.git
```

Before doing anything, install all the packages in the REQUIRE file. Then, run the following:

```
julia -e 'Pkg.clone(pwd()); Pkg.build(pwd()), Pkg.test(pwd(), coverage=true)'
```

This will make sure everything is built correctly.

# Using the package

Although more documentation is coming soon, this package supports alinging either one
pair or many pairs of sequences with the Smith-Waterman algorithm. There are also
functions to calculate ROC curves as well as optimize input score matrices. For now,
try running this command to align your favorite two sequences:

```
julia src/run.jl align -s scoring/BLOSUM50 -g 5 -e 3 sequences/test1.fa sequences/test2.fa output.txt
```

This will run the Smith Waterman algorithm to align "test1.fa" and "test2.fa" with a gap
penalty of 5 and extension penalty of 3. The scoring matrix can be specified with the
-s flag. The alignment and score are written out to `output.txt`.

You can also align multiple pairs of sequences and store the alignments in some file:

```
julia src/run.jl align_many -s scoring/BLOSUM50 -g 5 -e 3 scoring/Pospairs.txt output.txt
```

The package also supports functionality for outputting data for a ROC curve:

```
julia src/run.jl roc -s scoring/BLOSUM50 -g 5 -e 3 scoring/Pospairs.txt scoring/Negpairs.txt roc.output.txt
```

Finally, you can optimize a given substitution matrix:

```
julia src/utils/optimize_score_matrix scoring/MATIO .1 100 pos_alignments neg_alignments
```

for a simulated annealing parameter `T = .1` and for `N = 100` iterations. Pass static positive
alignments in `pos_alignments` and static negative alignments in `neg_alignments`.
