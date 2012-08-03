\name{getAE}
\alias{getAE}
\docType{data}
\title{ Download MAGE-TAB files from ArrayExpress in a specified directory }
\description{
  \code{getAE} downloads and extracts the MAGE-TAB files from an
  ArrayExpress dataset.
}
\usage{
getAE(input, path = getwd(), type = "raw", extract = TRUE)
}
\arguments{
  \item{input}{ is an ArrayExpress identifier. }
  \item{path}{ is the name of the directory in which the files
    downloaded on the ArrayExpress repository will be extracted.} 
  \item{type}{ can be 'raw' to download and extract only the raw data,
    'processed' to download and extract only the processed data or
    'full' to have both raw and processed data.}
  \item{extract}{ if FALSE, the files are not extracted from the zip archive.}
}
\value{
  \code{ A list with the names of the files that
    have been downloaded and extracted. }
}
\seealso{\code{\link[ArrayExpress]{ArrayExpress}},
  \code{\link[ArrayExpress]{ae2bioc}},
  \code{\link[ArrayExpress]{getcolproc}},
  \code{\link[ArrayExpress]{procset}}}
\author{
  Ibrahim Emam, Audrey Kauffmann 
  
  Maintainer: <iemam@ebi.ac.uk>  
}
\examples{
mexp1422 = getAE("E-MEXP-1422", type = "full")

## Build a an ExpressionSet from the raw data
MEXP1422raw = ae2bioc(mageFiles = mexp1422)

## Build a an ExpressionSet from the processed data
cnames = getcolproc(mexp1422)
MEXP1422proc = procset(mexp1422, cnames[2])

}
\keyword{datasets}
