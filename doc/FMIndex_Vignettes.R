## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(FMIndex)

## ----results="hide"-----------------------------------------------------------
#Specify the FASTA file, the parameters and the output directory
file <- system.file("extdata", "prova.fasta", package = "FMIndex")
step_tally <- 10
k <- 10
dir <- "./Prova_Results"

# Create the FM-Index in a specified directory
FM_index(file,10,10,dir)

## -----------------------------------------------------------------------------
# Creates the Suffix Array of a string
DNA_string <- "ACCTGGAC"
SA <- create_SA(DNA_string)
# To select a subset of the suffix array
# NOTE: this is already implemented in the FM_Index function
k <- 2 # set a value in line with the sequence length
SSA <- SA[SA %% k == 0] # Extract 

print(SA)
print(SSA)

## -----------------------------------------------------------------------------
# Creates the BWT of a DNA_string
DNA_string <- "ACCTGGAC"
bwt <- createBWT(DNA_string)
print(bwt)

## -----------------------------------------------------------------------------
DNA_string <- "ACCTGGAC"

#creates the sparse tally matrix of the BWT of the DNA string
step <- 2
sparse_tally <- sparse_tally(DNA_string,step)
print(sparse_tally)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

