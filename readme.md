# Advanced Methods in Bioinformatics - Tomás Guija Valiente & Benjamin Oberthür

## Content

- `./Report.pdf`: Review of the '**Comparing sequences without using
alignments: application to HIV/SIV subtyping** from *Didier et. al.*, 2007'. Source files from the report available in `./Report`
- `./data`:
  - `/Dissimilarity Matrices`: Dissimilarity matrices of the two sets of sequences
  - `/Newick trees`: Newick representation of the different clustering trees (from the original sequences and from bootstraps)
  - `/Phylip output`: Output files from the `Consensus.exe` program of `Phylip` (available in `./Resources/phylip-3.698`)
  - `/Rewritten sequences`: Rewritting of the sequences with the class identifiers
  - `/Tree visualization`: PNG files of the different trees (generated from the Newick representations with the following web app: http://www.trex.uqam.ca/index.php?action=newick)
- `./Resources`
  - `/LocDec - provided by Prof. Didier`: Sources files provided by Prof. Didier on our request, in which we found the sequence set we used in our experimentations
  - `/phylip-3.698`: Set of applications to deal with clustering trees from their Newick representation (retrieved from https://evolution.genetics.washington.edu/phylip.html)
  - `/Didier et al - 2007`: Article on which we base our review
  - `/German_Sequences.fasta`: Set of 50 sequences retrived from https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html
  - `/nef.fsa`: Copy of the sequence set from `/LocDec - provided by Prof. Didier`
- `./src`: Source files for the execution of our experimentation

> :warning: **WARNING**: The source file `./src/N_Local_Decoding.cpp` **works perfectly under Windows**, but seems to have **problems reading rewritten files under Linux Mint**.
>
> After multiple attemps on solving the problem, we decided to spend the remaining time on perfecting and completing our working code under the working environment.
>
> If any problem occurs under your computing environment, feel free to contact at: tomas.guija.valiente@ulb.be


