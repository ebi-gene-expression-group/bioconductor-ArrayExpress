\name{magetab2bioc}
\alias{magetab2bioc}
\docType{data}
\title{ Convert MAGE-TAB data into a Bioconductor object }
\description{
  \code{magetab2bioc} converts local MAGE-TAB files into a AffyBatch, an
  ExpressionSet or a NChannelSet.
}
\usage{
magetab2bioc(rawfiles, sdrf, idf, adf, path = ".", rawcol = NULL, save = TRUE)
}
\arguments{
  \item{rawfiles}{ all the expression files to use to create the object.
  The content of the raw.zip MAGE-TAB file.}
  \item{sdrf}{ the sdrf file from MAGE-TAB. }
  \item{idf}{ the idf file from MAGE-TAB. }
  \item{adf}{ the adf file from MAGE-TAB. }
  \item{path}{ is the name of the directory in which the files downloaded on the ArrayExpress repository will be extracted.}
  \item{rawcol}{ by default, the columns are automatically selected
    according to the scanner type. If the scanner is unknown or if the
    user wants to use different columns than the default, the argument
    'rawcol' can be set. For two colour arrays it must be a list with
    the fields 'R', 'G', 'Rb' and 'Gb' giving the column names to be
    used for red and green foreground and background. For one colour
    arrays, it must be a character string with the column name to be
    used. These column names must correspond to existing column names of
    the expression files.}
  \item{save}{ if TRUE, the files downloaded from the database will not be deleted from path after executing the function.}
}
\value{
  An object of class
  \code{\link[Biobase:class.AffyBatch]{AffyBatch}} or 
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} or 
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}} with the raw expression
  values in the assayData of the object, the information
  contained in the .sdrf file in the phenoData, the adf file content in
  the featureData and the idf file content in the experimentData.

  If several array designs are used in the data set, the output is a
  list with an object for each array design.  
}
\seealso{\code{\link{queryAE}}, \code{\link{getAE}}, \code{\link{ArrayExpress}}}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\examples{
MEXP1422 = getAE("E-MEXP-1422")
MEXP1422set = magetab2bioc(rawfiles = MEXP1422$rawfiles, sdrf =
  MEXP1422$sdrf, idf = MEXP1422$idf, adf = MEXP1422$adf)
}
\keyword{datasets}
