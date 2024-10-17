# FMIndex_package

For any FASTA file (containing a single DNA sequence) specified by the user, the R package will:
 1) Compute all necessary data structures of which the FM Index is composed
    - Suffix Array (SA)
    - Burrow-Wheeler Transform (BWT)
    - Character counts
    - Tally Matrix
 2) Write the individual FM data structures as individual .csv files into a user-specified folder. A compressed version of the tally matrix and the SA is saved.


See vignettes and manual for further information.


