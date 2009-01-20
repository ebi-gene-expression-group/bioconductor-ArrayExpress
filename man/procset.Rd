\name{procset}
\alias{procset}
\docType{data}
\title{ Convert processed MAGE-TAB files into Bioconductor object}
\description{
  \code{procset} converts local MAGE-TAB files into an ExpressionSet.
 }

\usage{
procset(files, procol, path, sdrf, adf, idf)
}

\arguments{
  \item{files}{ is the name of the processed MAGE-TAB file to be read. }
  \item{procol}{ the name of the column to be extracted from the
  file. Obtained using \code{\link{getcolproc}}.}
  \item{path}{ is the name of the directory where to find this file.}
  \item{sdrf}{ the sdrf file from MAGE-TAB. }
  \item{adf}{ the adf file from MAGE-TAB. }
  \item{idf}{ the idf file from MAGE-TAB. }
}

\seealso{\code{\link{ArrayExpress}}, \code{\link{getAE}},
  \code{\link{magetab2bioc}}, \code{\link{getcolproc}}, \code{\link{queryAE}}}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\examples{
files = getAE("E-MEXP-1422", type = "processed")
cnames = getcolproc(files$procfiles, path = ".")
MEXP1422 = procset(files$procfiles, cnames[2], ".", files$sdrf, files$adf, files$idf)
}
\keyword{datasets}
