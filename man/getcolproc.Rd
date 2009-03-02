\name{getcolproc}
\alias{getcolproc}
\docType{data}
\title{ Give the possible column names to from processed MAGE-TAB files }
\description{
  \code{getcolproc} extracts the column names from processed MAGE-TAB and
  return them, to be given as an input of the function \code{procset} 
}

\usage{
getcolproc(procfile, path = ".")
}

\arguments{
  \item{procfile}{ is the name of the processed MAGE-TAB file to be read. }
  \item{path}{ is the name of the directory where to find this file.}
}

\seealso{\code{\link{ArrayExpress}}, \code{\link{getAE}},
  \code{\link{magetab2bioc}}, \code{\link{procset}}, \code{\link{queryAE}}}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\examples{
mexp1422 = getAE("E-MEXP-1422", type = "processed")
cnames = getcolproc(mexp1422$procfile, path = ".")
}
\keyword{datasets}
