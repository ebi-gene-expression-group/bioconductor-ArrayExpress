\name{procset}
\alias{procset}
\docType{data}
\title{ Convert processed MAGE-TAB files into Bioconductor object}
\description{
  \code{procset} converts local MAGE-TAB files into an ExpressionSet.
 }

\usage{
procset(files, procol)
}

\arguments{
  \item{files}{ is the list with the names of the processed, the sdrf,
  the adf and the idf files and the path of the data. }
  \item{procol}{ the name of the column to be extracted from the
  file. Obtained using \code{\link{getcolproc}}.}
}

\seealso{\code{\link{ArrayExpress}}, \code{\link{getAE}},
  \code{\link{magetab2bioc}}, \code{\link{getcolproc}}, \code{\link{queryAE}}}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\examples{
mexp1422 = getAE("E-MEXP-1422", type = "processed")
cnames = getcolproc(mexp1422)
MEXP1422 = procset(mexp1422, cnames[2])
}
\keyword{datasets}
