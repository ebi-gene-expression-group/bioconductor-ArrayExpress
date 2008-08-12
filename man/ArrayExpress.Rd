\name{ArrayExpress}
\alias{ArrayExpress}
\docType{data}
\title{ R objects from ArrayExpress database }
\description{
  \code{ArrayExpress} produces an AffyBatch, an ExpressionSet or an
  NChannelSet from a raw data set from the ArrayExpress database.
  \code{ArrayExpress} needs an Internet connection. The object cannot be
  created if multiple chip types are used in the same data set.
}
\usage{
ArrayExpress(input, tempoutdir = ".", save = FALSE, columns = NULL)
}
\arguments{
  \item{input}{ an ArrayExpress identifier of an experiment for which
    the raw data are available. The package currently covers 96% of the
    raw data sets of the database.}
  \item{tempoutdir}{ the name of the directory in which the files
    downloaded on the ArrayExpress repository will be extracted. The
    default is the current directory.}
  \item{save}{ if TRUE, the expression files downloaded from the
    database will not be deleted from tempoutdir after executing the function.}
  \item{columns}{ by default, the columns are automatically selected
    according to the scanner type. If the scanner is unknown or if the
    user wants to use different columns than the default, the argument
    'columns' can be set. For two colour arrays it must be a list with
    the fields 'R', 'G', 'Rb' and 'Gb' giving the column names to be
    used for red and green foreground and background. For one colour
    arrays, it must be a character string with the column name to be
    used. These column names must correspond to existing column names of
    the expression files.}
  }
\value{
  An object of class
  \code{\link[Biobase:class.AffyBatch]{AffyBatch}} or 
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} or 
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}} with the expression
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
