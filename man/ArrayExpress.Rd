\name{ArrayExpress}
\alias{ArrayExpress}
\docType{data}
\title{ S4 class objects from ArrayExpress database }
\description{
  \code{ArrayExpress} produces an AffyBatch from an Affymetrix ArrayExpress dataset.
  \code{ArrayExpress} needs an Internet connection. The object cannot be
  created if multiple chip types are used in the same data set. The
  object cannot be created if expression data files are not correct
  binary .cel files.
}
\usage{
ArrayExpress(input, tempoutdir = ".")
}
\arguments{
  \item{input}{ is an ArrayExpress identifier of an Affymetrix experiment.}
  \item{tempoutdir}{ is the name of the directory in which the files
    downloaded on the ArrayExpress repository will be extracted. The
    default is the current directory.}
}
\value{
  An object of class
  \code{\link[Biobase:class.AffyBatch]{AffyBatch}} with the expression
  values in the assayData of the object and with the information
  contained in the .sdrf file in the phenoData.
}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\examples{
  \dontrun{
ETABM25.affybatch = ArrayExpress("E-TABM-25")
print(ETABM25.affybatch)
sampleNames(ETABM25.affybatch)
colnames(pData(ETABM25.affybatch))
}
}
\keyword{datasets}
