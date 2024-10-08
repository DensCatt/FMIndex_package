% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FMIndex_Script.R
\name{create_SA}
\alias{create_SA}
\title{\strong{Generate a Suffix Array for a DNA String}}
\usage{
create_SA(string)
}
\arguments{
\item{string}{A character string representing a DNA sequence.}
}
\value{
A numeric vector representing the suffix array, where each number is
the starting position of a suffix of the string in lexicographical order.
}
\description{
This function generates a suffix array (SA) for a given DNA string.
The SA stores the starting positions of the suffixes of the string in
lexicographical order.
This data structure is a component of the FM Index and is essential for a
more efficient creation of the BWT (\code{\link{createBWT}})
}
\examples{

DNA_string <- "ACCTGGAC"
SA <- create_SA(DNA_string) #
# To select a subset of the suffix array
 # NOTE: this is already implemented in the FM_Index function
k <- 2 # set a value in line with the sequence length
SSA <- SA[SA \%\% k == 0]

print(SA)
print(SSA)

}
