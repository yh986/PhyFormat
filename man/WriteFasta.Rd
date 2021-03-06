% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FastaPhylip.R
\name{WriteFasta}
\alias{WriteFasta}
\title{FASTA file writer}
\usage{
WriteFasta(x, filename, multiline = 60)
}
\arguments{
\item{x}{The physeq object to read from}

\item{filename}{The name of the file to create}

\item{multiline}{How many characters per line, for an interleaved FASTA format. Use 'multiline = -1' for sequential format (no character limit). Default = 60}
}
\description{
Writes a fasta file from a physeq object.
}
\examples{
x <- ReadPhylip("MSA.phy") #Creates the physeq object 'x'
WriteFasta(x, "MSA.fa")
}
