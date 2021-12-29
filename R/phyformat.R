#' phyformat: A package to help bioinformatics data formatting
#'
#' phyformat provides the ability to read and write sequence data in FASTA and PHYLIP formats, including the ability to interconvert between the two.
#' It also includes a function to take a subset of a Newick-format phylogenetic tree.
#' 
#' @examples 
#' x <- ReadFasta("sequences.fa") #Creates a physeq object
#' x$Sequences[[1]] <- "PQHLIRV" #Edit a sequence in the physeq object
#' WriteFasta(x, "sequences2.fa") #Write to a fasta file
#' WritePhylip(x, "sequences2.phy") #Write to a PHYLIP file
#' @docType package
#' @name phyformat
NULL
#> NULL
