\name{getAE}
\alias{getAE}
\docType{data}
\title{ Download MAGE-TAB from ArrayExpress in a specified directory }
\description{
  \code{getAE} downloads and extracts the MAGE-TAB files from an
  ArrayExpress data set.
}
\usage{
getAE(input, path = ".", save = TRUE, type = "full", extract = TRUE)
}
\arguments{
  \item{input}{ is an ArrayExpress identifier. }
  \item{path}{ is the name of the directory in which the files downloaded on the ArrayExpress repository will be extracted.}
  \item{save}{ if TRUE, the files downloaded from the database will not be deleted from path after executing the function.}
  \item{type}{ can be 'raw' to download and extract only the raw data,
    'processed' to download and extract only the processed data or
    'full' to have both raw and processed data.}
  \item{extract}{ if FALSE, the files are not extracted from the zip archive.}
}
\value{
  \code{ A list with the names of the files that
    have been downloaded and extracted.}
}
\seealso{\code{\link{ArrayExpress}}, \code{\link{magetab2bioc}}, \code{\link{queryAE}}}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\examples{
TABM25 = getAE("E-TABM-25", type = "raw")
}
\keyword{datasets}
