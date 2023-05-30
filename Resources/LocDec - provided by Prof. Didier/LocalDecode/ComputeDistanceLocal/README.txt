usage: local [options] <input file> <output file>

The input file has to be in Fasta format.

Options:

-h                  output help

-f <f, s, m, x>     select the format of output for distances matrix
                    (r: raw, t: text table, p: phylip, n: nexus)

-s <d, r, p>        select the type of alphabet (d: DNA, r: RNA, p: Protein)
                    maybe useful to handle ambiguity characters

-d <l, s>           select the type of decoding:
                     * l: use local decoding
                     * s: use sequence of sliding blocks

-m <number>         select the type of dissimilarity:
                     * 0: trivial dist
                     * 1: Pham
-o <number>         set the order of local decoding/length of the blocks