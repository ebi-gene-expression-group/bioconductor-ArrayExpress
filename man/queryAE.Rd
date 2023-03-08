\name{queryAE}
\alias{queryAE}
\docType{data}
\title{ Query of the ArrayExpress collection in Biostudies }
\description{
  \code{queryAE} queries the ArrayExpress collection with keywords and
  give a dataframe with ArrayExpress identifiers and related
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
         such as the number of files, the release date on the database, 
         the title, the author and content.}
}
\seealso{\code{\link[ArrayExpress]{ArrayExpress}},
  \code{\link[ArrayExpress]{getAE}}} 
\author{
  Ibrahim Emam, Audrey Kauffmann
  
  Maintainer: Jose Marugan <jcmarca@ebi.ac.uk>  
}
\examples{
## To retrieve all the identifiers of pneumonia data sets
pneumo = queryAE(keywords = "pneumonia")

## To retrieve all the identifiers of pneumonia data sets studied in human
pneumoHS = queryAE(keywords = "pneumonia", species = "homo+sapiens")
}
\keyword{datasets}
