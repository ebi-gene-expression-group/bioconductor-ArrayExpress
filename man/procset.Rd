\name{procset}
\alias{procset}
\docType{data}
\title{ Convert processed MAGE-TAB files into a Bioconductor object}
\description{
  \code{procset} converts local MAGE-TAB files into an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.
 }

\usage{
procset(files, procol)
}

\arguments{
  \item{files}{ is the list with the names of the processed, the sdrf,
    the adf and the idf files and the path of the data as given by
    \code{\link{getAE}}. }
  \item{procol}{ the name of the column to be extracted from the
  file. Obtained using \code{\link{getcolproc}}.}
}

\seealso{\code{\link{queryAE}}, \code{\link{getAE}}, \code{\link{getcolproc}}}
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