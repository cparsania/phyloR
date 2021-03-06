% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blast_result_parser.R
\name{format_fasta_headers}
\alias{format_fasta_headers}
\title{Format fasta headers}
\usage{
format_fasta_headers(fasta_file = NULL, keep_alignemnt_coord = TRUE)
}
\arguments{
\item{fasta_file}{string denoting full path of a fasta file.}

\item{keep_alignemnt_coord}{logical (default : TRUE) decides whether alignment coordinates to keep in headers or not. When TRUE header must contains alignment coordinates in this format \code{:([:digit:]+-[:digit:]+)}. (e.g. CEJ90625.1:1-252)}
}
\value{
sequences as an object of class \code{Biostrings}
}
\description{
format the sequence headers for fasta file obtained from NCBI blast output
}
\examples{
\dontrun{
 f <- system.file("extdata","blast_output_01.fasta" ,package = "phyloR")
 ## existing headers
 Biostrings::readAAStringSet(f) \%>\% names() \%>\% head()
 ## new  headers
 f \%>\% format_fasta_headers() \%>\% names() \%>\% head()
}
}
