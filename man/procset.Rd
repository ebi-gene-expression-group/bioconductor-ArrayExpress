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
    \code{\link[ArrayExpress]{getAE}}. }
  \item{procol}{ the name of the column to be extracted from the
  file. Obtained using \code{\link[ArrayExpress]{getcolproc}}.}
}

\seealso{\code{\link[ArrayExpress]{queryAE}},
  \code{\link[ArrayExpress]{getAE}},
  \code{\link[ArrayExpress]{getcolproc}}} 
\author{
  Ibrahim Emam, Audrey Kauffmann 
  
  Maintainer: <iemam@ebi.ac.uk>  
}
\examples{
# An example can be found in the help of the getAE function.
}
\keyword{datasets}
