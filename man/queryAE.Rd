\name{queryAE}
\alias{queryAE}
\docType{data}
\title{ XML query of the ArrayExpress repository }
\description{
  \code{queryAE} queries the ArrayExpress database with keywords and
  give a dataframe with  ArrayExpress identifiers and related
  information, as an output.
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
  \code{ A dataframe with all the ArrayExpress dataset
         identifiers which correspond to the query in the first column.
         The following columns contain information about these datasets, 
         such as the number of raw files, the number of data processed, 
         the release date on the database, the pubmed ID, the species, 
         the experiment design and the experimental factors.}
}
\seealso{\code{\link{ArrayExpress}}, \code{\link{getAE}}}
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
