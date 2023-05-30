usage: compare [options] <input file1> <input file2>

returns the name of input file name followed by the Pearson coefficient
between the two distance (Mantel score)

The input files have to be in Nexus format.

Options:

-h                  output help

-o <output file>    write a table with in each line values for dist1 and dist2

-f <output file>    write distances matrix "dist1-dist2"

-f <f, s, m, x>     select the format of output for distances matrix
                    (r: raw, t: text table, p: phylip, n: nexus)
