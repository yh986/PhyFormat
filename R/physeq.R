#' physeq constructor
#' 
#' Creates a new physeq object (a list of sequences) from a vector of headings and sequences, which can then be saved as a FASTA or PHYLIP file.
#' 
#' @param headings Vector of sequence names
#' @param sequences Vector of nucleotide/amino acid sequences
#' @export
#' @examples 
#' x <- physeq(c("Seq 1", "Seq 2"), c("RRARRPR", "PKKKRKV"))
#' print(x)
#' x$Headings #Get or modify headings
#' x$Sequences #Get or modify sequences
#' WriteFasta(x, "Seqs.fa")
physeq <- function (headings, sequences) {
  return (validate_physeq(new_physeq(headings, sequences)))
}

#Internal constructor
#Just makes a list of the two input vectors and puts a class on it (S3)
new_physeq <- function (headings, sequences) {
  stopifnot(is.character(headings) || is.character(sequences))
  structure(list(Headings = headings, Sequences = sequences), class = "physeq")
}

#Validates a physeq to see that it makes sense
validate_physeq <- function (x) {
  unclass <- unclass(x)
  #Verify that there are as many headings as sequences
  
  if (length(unclass$Headings) != length(unclass$Sequences)) {
    stop("There must be as many headings as there are sequences")
  }
  #Additional criteria:
  #-Headings should be unique (unchecked)
  #-Sequences should have valid characters, depending on if it is amino acids or nucleotides (unchecked)
  
  x
}

#S3 method for displaying a physeq
#' @export
print.physeq <- function(x, ...) {
  numsequences <- length(x$Headings)
  cat(numsequences, "sequences. \n")
  
  printtruncate <- function(str, n) {
    if (nchar(str) > n) {
      paste0(substr(str, 1, n - 3), "...")
    }
    else
    {
      substr(str, 1, n)
    }
  }
  
  for (i in 1:min(numsequences, 5)) {
    cat(paste0("\"", printtruncate(x$Headings[i], 15), "\""), ":", printtruncate(x$Sequences[i], 20), "\n")
  }
  if (numsequences > 5) {
    cat("...")
  }
  
  invisible(x)
}
