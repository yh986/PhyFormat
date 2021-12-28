#' FASTA file reader
#'
#' Creates a new physeq object from a FASTA-format file.
#' @param filename The name of the file to be read
#' @export
#' @examples 
#' x <- ReadFasta("Seqs.fa")
#' x$Sequences #Get or modify sequences
ReadFasta <- function(filename) {
  
  if (!file.exists(filename)) {
    stop("Could not find the specified file")
  }
  
  fa <- readLines(filename)
  faheaderlines <- which(grepl("^>.*", fa))
  numberofseqs <- length(faheaderlines)
  separators <- c(faheaderlines, length(fa) + 1) #Get the line numbers corresponding to the beginnings and ends of sequences
  seqvector <- NULL
  headervector <- NULL
  for (i in 1:numberofseqs) {
    beginseq <- separators[i] + 1
    endseq <- separators[i + 1] - 1
    seq_i <- paste0(fa[beginseq:endseq], collapse = "")
    header_i <- substr(fa[separators[i]], 2, nchar(fa[separators[i]]))
    seqvector <- c(seqvector, seq_i)
    headervector <- c(headervector, header_i)
  }
  
  out_physeq <- new_physeq(headervector, seqvector)
  validate_physeq(out_physeq)
  return (out_physeq)
}

#' PHYLIP file reader
#' 
#' Creates a new physeq object from an existing PHYILP-format file.
#' @param filename The name of the file to be read
#' @export
#' @examples 
#' x <- ReadPhylip("MSA.phy")
#' x$Sequences #Get or modify sequences
ReadPhylip <- function(filename) {
  
  if (!file.exists(filename)) {
    stop("Could not find the specified file")
  }
  
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
  headers <- substr(rawheads, 1, n)
  headers <- trimws(headers, which = "right")
  
  out_physeq <- new_physeq(headers, sequences)
  validate_physeq(out_physeq)
  return (out_physeq)
}


#' FASTA file writer
#' 
#' Writes a fasta file from a physeq object.
#' @param x The physeq object to read from
#' @param filename The name of the file to create
#' @param multiline How many characters per line, for an interleaved FASTA format. Use 'multiline = -1' for sequential format (no character limit). Default = 60
#' @export
#' @examples 
#' x <- ReadPhylip("MSA.phy") #Creates the physeq object 'x'
#' WriteFasta(x, "MSA.fa")
WriteFasta <- function(x, filename, multiline = 60) {
  x <- unclass(x)
  headers <- x[["Headers"]]
  sequences <- x[["Sequences"]]
  
  writestring <- NULL
  for (i in 1:length(headers)) {
    writestring <- c(writestring, paste0(">", headers[i]))
    
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

#' PHYLIP-format file sequence alignment writer
#' 
#' Writes a PHYLIP multiple sequence alignment file from a physeq object. If sequences are not the same length, there will be an error
#' @param x The physeq object to read from
#' @param filename The name of the file to create
#' @param blocks The number of sequences per block. Default = 5
#' @param blocksize The number of characters per sequence per block. Default = 10
#' @examples 
#' x <- ReadFasta("MSA.fa") #Creates the physeq object 'x'
#' WritePhylip(x, "MSA.phy")
#' y <- physeq(c("Seq 1", "Seq 2"), c("RRARRPR", "PKKKRKV"))
#' WritePhylip(y, "Seqs.phy")
#' @export
WritePhylip <- function(x, filename, blocks = 5, blocksize = 10) {
  x <- unclass(x)
  headers <- x[["Headers"]]
  sequences <- x[["Sequences"]]
  
  if (!all(nchar(headers) == nchar(headers[1]))) {
    stop("Sequences are not the same length")
  }
  
  #Write the first line
  N <- length(headers)
  len <- nchar(sequences[[1]])
  writestring <- paste0(N, " ", len)
  #Write the first paragraph
  paragraph <- paste0(headers, "          ")
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
