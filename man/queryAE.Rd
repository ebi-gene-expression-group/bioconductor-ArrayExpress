\name{queryAE}
\alias{queryAE}
\docType{data}
\title{ XML query of the ArrayExpress repository }
\description{
  \code{queryAE} queries ArrayExpress with keywords and give a list of
  ArrayExpress identifiers and their availability as raw and/or
  processed data, as an output.
}
\usage{
queryAE(keywords = NULL, species = NULL)
}
\arguments{
  \item{keywords}{ the keyword(s) of interest. To use several words,
    they must be separated by a "+" as shown in the examples. }
  \item{species}{ the specie(s) of interest.}
}
\value{
  \code{ A dataframe with character strings with all the ArrayExpress dataset
    identifiers which correspond to the query in the first column. The
    two following columns are filled with 'Yes' and 'No' according to
    the existence of raw and processed data.}
}
\seealso{\code{\link{ArrayExpress}}, \code{\link{getAE}}, \code{\link{magetab2bioc}}}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\examples{
## To retrieve all the identifiers of cancer data sets
Cancer = queryAE(keywords = "cancer")

## To retrieve all the identifiers of cancer data sets studied in human
CancerHS = queryAE(keywords = "cancer", species = "homo+sapiens")

## To retrieve all the identifiers of breast cancer data sets studied in human
BreastCancerHS = queryAE(keywords = "breast+cancer", species = "homo+sapiens")
}
\keyword{datasets}
