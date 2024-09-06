#' \strong{Generate a Suffix Array for a DNA String}
#'
#' This function generates a suffix array (SA) for a given DNA string.
#' The SA stores the starting positions of the suffixes of the string in
#' lexicographical order.
#' This data structure is a component of the FM Index and is essential for a
#' more efficient creation of the BWT (\code{\link{createBWT}})
#'
#' @param string A character string representing a DNA sequence.
#' @return A numeric vector representing the suffix array, where each number is
#' the starting position of a suffix of the string in lexicographical order.
#'
#' @examples
#'
#' DNA_string <- "ACCTGGAC"
#' SA <- create_SA(DNA_string) #
#' # To select a subset of the suffix array
#'  # NOTE: this is already implemented in the FM_Index function
#' k <- 2 # set a value in line with the sequence length
#' SSA <- SA[SA %% k == 0]
#'
#' print(SA)
#' print(SSA)
#'
#' @export
create_SA <- function(string) {

    string <- paste0(string, "$")
    n <- nchar(string)

    suffix_array <- seq_len(n) #Generates the indexes of the suffixes
    suffixes <- vapply(suffix_array, function(i) substring(string, i, n),
                                                            character(1))
    sorted_indices <- order(suffixes)
    suffix_array <- suffix_array[sorted_indices]
    return(suffix_array)
}


#' \strong{Creates a Burrows-Wheeler Transform of a DNA string}
#'
#' This function takes a DNA string as input, containing the bases A,C,G and T.
#' Then creates a Burrows-Wheeler Transformation (BWT) of that string through
#' the usage of a suffix array.
#' The BWT is the last column of the Burrows-Wheeler Matrix (BWM).
#' It gives in output a vector of 5 character containing the BWT
#' (as the last element) and 4 strings representing a compressed version of the
#' first column of the Burrows_Wheeler Matrix (BWM).
#' These strings correspond to the counts of the bases A, C, G, and T
#' (in that order) in the first column of the BWM, respectively.
#'
#' @param string a character string representing a DNA sequence
#' @return a vector of 5 strings.
#'        The first 4 elements represent the count of A's, C's, G's, and T's.
#'        The last element is the BWT of the DNA string.
#' @seealso \code{\link{create_SA}} for generating the suffix array used in this
#'        function.
#' @examples
#' DNA_string <- "ACCTGGAC"
#' bwt <- createBWT(DNA_string)
#' print(bwt)
#'
#' @export
createBWT <- function(string) {

    suffix_array <- create_SA(string)

    string <- paste0(string, "$")
    n <- nchar(string)

    bwt <- vapply(suffix_array, function(i) {
        if (i == 1) {
        # if the suffix is the entire string, take the last char
            return(substring(string, n, n))}
        else {
        # else take the previous one
            return(substring(string, i - 1, i - 1))
        }
        },character(1))

    first_column <- vapply(suffix_array, function(i) {
        return(substring(string, i, i))
        }, character(1))

    first_col_counts <- as.numeric(table(factor(first_column,
                        levels = c("A", "C", "G", "T"))))
    bwt <- paste(bwt, collapse = "")

    return(c(first_col_counts,BWT=bwt))
}


#' \strong{Generate a Sparse Tally Matrix of a Burrows-Wheeler Transform}
#'
#' This function generates a sparse tally matrix of the Burrows-Wheeler
#' Transform (BWT) of a given DNA string.
#' The BWT is generated inside the function, so the input is just the DNA string
#' and a step.
#' The sparse tally matrix is useful for efficiently storing cumulative
#' character counts at specific intervals, specified by the step parameter,
#' which can be used for further applications, such as FM-index creation.
#'
#' @param string A character string representing a DNA sequence, containing the
#'        bases A, C, G, and T.
#'        The BWT is generated within the function.
#' @param step An integer specifying the step size used for selecting rows from
#'        the tally matrix to create its sparse version containing just the
#'        selected rows
#' @return A sparse tally matrix representing cumulative counts of each base
#'        (A, C, G, T) in the BWT at the specified intervals.
#'        The row names correspond to the positions in the original BWT.
#' @details The function first computes the BWT of the input string
#'          using \link{createBWT}. It then calculates cumulative counts for
#'          each unique character in the BWT.
#'          The resulting sparse matrix stores these counts only at the
#'          intervals defined by `step`, making it more memory-efficient while
#'          retaining essential information for sequence analysis.
#' @examples
#' DNA_string <- "ACGCGGATTGTGATC"
#' sparse_matrix <- sparse_tally(DNA_string, 2)
#' print(sparse_matrix)
#'
#' @seealso \code{\link{createBWT}} for generating the Burrows-Wheeler Transform
#'
#' @export
sparse_tally <- function(string, step) {

    chars <- strsplit(createBWT(string)["BWT"], NULL)[[1]]

    unique_chars <- unique(chars)

    presence_matrix <- vapply(unique_chars, function(char)
                            as.numeric(chars == char),
                            FUN.VALUE = numeric(length(chars)))

    result_matrix <- apply(presence_matrix, 2, cumsum)

    selected_rows <- seq(from = 1, to = nrow(result_matrix), by = step)
    sparse_matrix <- result_matrix[selected_rows, , drop = FALSE]

    rownames(sparse_matrix) <- selected_rows
    colnames(sparse_matrix) <- unique_chars

    return(sparse_matrix)
}



#' \strong{Create a compressed FM-index from a DNA Sequence}
#'
#' This function creates an FM-index for a given DNA sequence in FASTA format.
#' The FM-index is a data structure that is useful for an efficient pattern
#' search and storage. It consists of the \emph{Burrows-Wheeler Transform (BWT),
#' a sparse tally matrix, and a subset suffix array (SSA).}
#' The resulting data structures are saved as CSV files
#' in a specified directory.
#' This compressed version of the FM-Index can potentially store an entire human
#' genome in less than 2 GB
#' Increasing \code{step_tally} and \code{k}, leads to a more efficient storage,
#' on the other hand, it could lead to longer search times for patterns.
#'
#' @param fasta_file A string representing the path to a FASTA file containing
#'        a DNA sequence.
#' @param step_tally An integer specifying the step size for generating the
#'        sparse tally matrix.
#'        Default is 32.
#' @param k An integer specifying the step size for selecting elements from the
#'        suffix array to create the SSA.
#'        Default is 32.
#' @param output A string specifying the folder that will contain the output
#'        files (BWT, sparse tally, SSA) should be saved.
#'        If it does not exist, will be created.
#'        Default is ./FM_index_Results
#' @return A message of success and the location of the saved FM-Index files.
#'         If an error occurs, an error message is returned, and the function
#'         execution is stopped.
#'
#' @import GenomicRanges
#' @importFrom Biostrings readDNAStringSet
#'
#' @examples
#'
#' file <- system.file("extdata", "prova.fasta", package = "FMIndex")
#' FM_index(file, step_tally = 10, k = 10)
#'
#' @seealso \code{\link{createBWT}}, \code{\link{sparse_tally}},
#'           \code{\link{create_SA}}
#' @export
FM_index <- function(fasta_file,step_tally=32,k=32,output="./FM_index_Results"){
    if (!endsWith(fasta_file,".fasta")) {
        return("ERROR: input DNA sequence must be in FASTA format")
    }
    if (!file.exists(fasta_file)) {
        return(paste0("ERROR: cannot find file '",fasta_file,"'"))
    }

    DNA_seq <- Biostrings::readDNAStringSet(fasta_file)

    if (step_tally > Biostrings::width(DNA_seq)/2) {
        return(paste0("ERROR: step_tally value (",step_tally,")
                        greater than half the DNA sequence length"))
    }
    if (k > Biostrings::width(DNA_seq)/2) {
        return(paste0("ERROR: k value (",k,")
                        greater than half the DNA sequence length"))
    }
    if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
    else {
        return(paste("ERROR: output folder",output,"already exists. ",
                        "Choose a different directory/folder to
                        avoid overwriting"))}

    string <- as.character(DNA_seq)[[1]]

    bwt <- createBWT(string)
    tally <- sparse_tally(string,step_tally)
    SA <- create_SA(string)

    SSA <- SA[SA %% k == 0] #reduce the memory required

    write.table(bwt, file = file.path(output, "BWT.csv"),
                row.names = c("A:","C:","G:","T:","BWT:"),
                quote=FALSE,col.names=FALSE,sep=";")
    write.table(tally, file = file.path(output, "Sparse_Tally.csv"),
                row.names = FALSE,sep=";",quote=FALSE)
    write.table(SSA, file = file.path(output, "Subset_Suffix_Array.csv"),
                row.names=FALSE,col.names = FALSE,sep=";")

    return(paste0("FM-Index of the DNA sequence '",fasta_file,"'created.
                The compressed data structures are saved in",output))}

