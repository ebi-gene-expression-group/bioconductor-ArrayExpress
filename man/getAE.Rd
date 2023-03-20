\name{getAE}
\alias{getAE}
\docType{data}
\title{ Download MAGE-TAB files from the ArrayExpress collection in a specified directory }
\description{
  \code{getAE} downloads and extracts the MAGE-TAB files from an
  ArrayExpress dataset.
}
\usage{
getAE(accession, path = getwd(), type = "full", extract = TRUE, sourcedir = path, overwrite = FALSE)
}
\arguments{
  \item{accession}{ is an ArrayExpress experiment identifier. }
  \item{path}{ is the name of the directory in which the files
    downloaded on the ArrayExpress repository will be extracted.} 
  \item{type}{ can be 'raw' to download and extract only the raw data,
    'processed' to download and extract only the processed data or
    'full' to have both raw and processed data.}
  \item{extract}{ if FALSE, the files are not extracted from the zip archive.}
  \item{sourcedir}{ when local = TRUE, files will be read from this directory.}
  \item{overwrite}{ if TRUE, overwrite files if they already exist in path, default FALSE.}
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
  
  Maintainer: Jose Marugan <jcmarca@ebi.ac.uk>  
}
\examples{
mexp21 = getAE("E-MEXP-21", type = "full")

## Build a an ExpressionSet from the raw data
MEXP21raw = ae2bioc(mageFiles = mexp21)

## Build a an ExpressionSet from the processed data
cnames = getcolproc(mexp21)
MEXP21proc = procset(mexp21, cnames[2])

}
\keyword{datasets}
