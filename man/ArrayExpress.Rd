\name{ArrayExpress}
\alias{ArrayExpress}
\docType{data}
\title{ S4 class objects from ArrayExpress database }
\description{
  \code{ArrayExpress} produces an AffyBatch from an Affymetrix ArrayExpress dataset.
  \code{ArrayExpress} needs Internet connection.
}
\usage{
ArrayExpress(input, tempoutdir)
}
\arguments{
  \item{input}{ is an ArrayExpress identifier.}
  \item{tempoutdir}{ is the name of the directory in which the files downloaded on the ArrayExpress repository will be extracted.}
}
\value{
  A S4 class object of class
  \code{\link[Biobase:class.AffyBatch]{AffyBatch}} if the input data set
  is an Affymetrix experiment,
}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\examples{
  \dontrun{
EAFMX4.affybatch = ArrayExpress("E-AFMX-4", tempoutdir = ".")
}
}
\keyword{datasets}
