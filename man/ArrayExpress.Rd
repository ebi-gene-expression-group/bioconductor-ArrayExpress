\name{ArrayExpress}
\alias{ArrayExpress}
\docType{data}
\title{ R objects from ArrayExpress database }
\description{
  \code{ArrayExpress} produces an AffyBatch, an ExpressionSet or an
  NChannelSet from a raw data set from the ArrayExpress
  database. \code{ArrayExpress} needs an Internet connection.
}
\usage{
ArrayExpress(input, path = ".", save = TRUE, rawcol = NULL)
}
\arguments{
  \item{input}{ an ArrayExpress identifier.
    The package currently covers 96\% of the raw data sets of the
    database and is not currently dealing with the processed data.}
  \item{path}{ the name of the directory in which the files
    downloaded on the ArrayExpress repository will be extracted. The
    default is the current directory.}
  \item{save}{ if TRUE, the files downloaded from the
    database will not be deleted from path after executing the function.}
  \item{rawcol}{ by default, for the raw data, the columns are
    automatically selected
    according to the scanner type. If the scanner is unknown or if the
    user wants to use different columns than the default, the argument
    'rawcol' can be set. For two colour arrays it must be a list with
    the fields 'R', 'G', 'Rb' and 'Gb' giving the column names to be
    used for red and green foreground and background. For one colour
    arrays, it must be a character string with the column name to be
    used. These column names must correspond to existing column names of
    the expression files.}
  }
\value{
  The output is an object of class
  \code{\link[Biobase:class.AffyBatch]{AffyBatch}} or 
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} or 
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}} with the raw expression
  values in the assayData of the object, the information
  contained in the .sdrf file in the phenoData, and the idf file content
  in the experimentData.

  If several array designs are used in the data set, the output is a
  list with an object for each array design.
}
\seealso{\code{\link{queryAE}}, \code{\link{getAE}}, \code{\link{magetab2bioc}}}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\examples{
ETABM25.affybatch = ArrayExpress(input = "E-TABM-25")
print(ETABM25.affybatch)
sampleNames(ETABM25.affybatch)
colnames(pData(ETABM25.affybatch))
}
\keyword{datasets}
