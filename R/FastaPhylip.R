#' @export
new_phytools_msa <- function (headings, sequences) {
  stopifnot(is.character(headings) || is.character(sequences))
  structure(list(Headings = headings, Sequences = sequences), class = "phytools_msa")
}

#' @export
validate_phytools_msa <- function (x) {
  unclass <- unclass(x)
  #Verify that there are as many headings as sequences
  
  if (length(unclass$Headings) != length(unclass$Sequences)) {
    stop("There must be as many headings as there are sequences")
  }
  #Additional criteria:
  #-Headings should be unique (unchecked)
  #-Sequences should have valid characters (unchecked)
  #-Sequences should have the same length
  
  sequences <- unclass$Sequences
  samelengths <- all(nchar(sequences) == nchar(sequences[1]))
  
  if (is.na(samelengths)) {
    stop("Undefined sequence")
  }
  
  if (!isTRUE(samelengths)) {
    stop("Alignment lengths have different lengths")
  }
  
  x
}

print.phytools_msa <- function(x, ...) {
  numsequences <- length(x$Headings)
  cat("MSA with", numsequences, "sequences. \n")
  
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
    cat(paste0("\"", printtruncate(x$Headings[i], 15), "\""), ":", printtruncate(x$Sequences[i], 20))
  }
  if (numsequences > 5) {
    cat("...")
  }
  
  invisible(x)
}

#' @export
ReadFasta <- function(filename) {
  fa <- readLines(filename)
  faheaderlines <- which(grepl("^>.*", fa))
  numberofseqs <- length(faheaderlines)
  separators <- c(faheaderlines, length(fa) + 1)
  seqvector <- NULL
  headingvector <- NULL
  for (i in 1:numberofseqs) {
    beginseq <- separators[i] + 1
    endseq <- separators[i + 1] - 1
    seq_i <- paste0(fa[beginseq:endseq], collapse = "")
    heading_i <- substr(fa[separators[i]], 2, nchar(fa[separators[i]]))
    seqvector <- c(seqvector, seq_i)
    headingvector <- c(headingvector, heading_i)
  }
  
  out_phytools_msa <- new_phytools_msa(headingvector, seqvector)
  validate_phytools_msa(out_phytools_msa)
  return (out_phytools_msa)
}

#' @export
ReadPhylip <- function(filename) {
  phylip <- readLines(filename)
  
  N <- as.numeric(sub("^[^0-9]*([0-9]+)[^0-9]*([0-9]+)[^0-9]*$", "\\1", phylip[1]))
  len <- as.numeric(sub("^[^0-9]*([0-9]+)[^0-9]*([0-9]+)[^0-9]*$", "\\2", phylip[1]))
  seqlines <- phylip[grepl("[-A-Za-z*]", phylip)]
  #Make N concatenated sequences
  concatlines <- rep("", N)
  for (i in 1:(length(seqlines)/N)) {
    beginparagraph <- ((i - 1) * N) + 1
    endparagraph <- i * N
    concatlines <- paste0(concatlines, seqlines[beginparagraph:endparagraph])
  }
  
  rawheads <- substr(concatlines, 1, 10)
  tails <- substr(concatlines, 11, nchar(concatlines[1]))
  tails_nospace <- gsub("[^-A-Za-z*]", "", tails)
  
  #Determine how much of the headers' ending to include in the actual sequence
  lenaccounted <- nchar(tails_nospace[1])
  n <- 10 #Last position of headers
  while (lenaccounted < len) {
    char_n <- substr(rawheads, n, n)
    if (grepl("[-A-Za-z*]", char_n)) {
      lenaccounted <- lenaccounted + 1
    }
    
    if (n <= 1) {
      break
    }
    else
    {
      n <- n - 1
    }
  }
  sequences <- paste0(substr(rawheads, n + 1, 10), 
                      tails_nospace)
  headings <- substr(rawheads, 1, n)
  headings <- trimws(headings, which = "right")
  
  out_phytools_msa <- new_phytools_msa(headings, sequences)
  validate_phytools_msa(out_phytools_msa)
  return (out_phytools_msa)
}

#E.g. multiline = 60 to break up sequences every 60 positions
#multiline = -1 for no multiline
#' @export
WriteFasta <- function(x, filename, multiline = 60) {
  x <- unclass(x)
  headings <- x[["Headings"]]
  sequences <- x[["Sequences"]]
  
  writestring <- NULL
  for (i in 1:length(headings)) {
    writestring <- c(writestring, paste0(">", headings[i]))
    
    if (multiline <= 0)
    {
      writestring <- c(writestring, sequences[i])
    }
    else
    {
      n <- 0
      while (TRUE) {
        substringbegin <- n * multiline + 1
        substringend <- (n + 1) * multiline
        nextsubstring <- substr(sequences[i], substringbegin, substringend)
        writestring <- c(writestring, nextsubstring)
        
        if (substringend < nchar(sequences[i]))
        {
          #Add new line
          n <- n + 1
        }
        else
        {
          break
        }
      }
    }
  }
  fileConn <- file(filename)
  writeLines(writestring, fileConn)
  close(fileConn)
}

#' @export
WritePhylip <- function(x, filename, blocks = 5, blocksize = 10) {
  x <- unclass(x)
  headings <- x[["Headings"]]
  sequences <- x[["Sequences"]]
  #Write the first line
  N <- length(headings)
  len <- nchar(sequences[[1]])
  writestring <- paste0(N, " ", len)
  #Write the first paragraph
  paragraph <- paste0(headings, "          ")
  paragraph <- substr(paragraph, 1, 10)
  paragraph <- paste0(paragraph, " ")
  n <- 0 #Block counter
  substringbegin <- 1
  substringend <- blocksize
  
  while(TRUE) {
    paragraph <- paste0(paragraph, 
                        substr(sequences, substringbegin, substringend)
    )
    if (n < (blocks - 1) && substringend < len)
    {
      #Write the next block
      paragraph <- paste0(paragraph, " ")
      n <- n + 1
      substringbegin <- substringbegin + blocksize
      substringend <- substringend + blocksize
    }
    else
    {
      break
    }
  }
  
  writestring <- c(writestring, paragraph) #Add paragraph to growing string
  
  #Write the next paragraphs
  while (TRUE) {
    #Is next paragraph needed?
    if (substringend < len) {
      writestring <- c(writestring, "") #Separator for next paragraph
      paragraph <- rep("", N)
      n <- 0
      substringbegin <- substringbegin + blocksize
      substringend <- substringend + blocksize
      while (TRUE) {
        paragraph <- paste0(paragraph, 
                            substr(sequences, substringbegin, substringend)
        )
        
        
        if (n < (blocks - 1) && substringend < len)
        { 
          #Add block
          paragraph <- paste0(paragraph, " ")
          n <- n + 1
          substringbegin <- substringbegin + blocksize
          substringend <- substringend + blocksize
          #Loop
        }
        else
        {
          break
        }
        
      }
      #After inner while loop is done:
      writestring <- c(writestring, paragraph) #Add paragraph to growing string
    }
    else
    {
      break
    }
  }
  fileConn <- file(filename)
  writeLines(writestring, fileConn)
  close(fileConn)
}
