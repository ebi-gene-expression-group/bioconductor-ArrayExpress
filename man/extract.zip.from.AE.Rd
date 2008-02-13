\name{extract.zip.from.AE}
\alias{extract.zip.from.AE}
\docType{data}
\title{ Unzip archives in a specified directory }
\description{
  \code{extract.zip} extracts the files from a .zip archive in
  ArrayExpress repository.
}
\usage{
extract.zip.from.AE(input, unzip = getOption("unzip"), tempoutdir = tempdir())
}
\arguments{
  \item{input}{ An ArrayExpress identifier. }
  \item{unzip}{ character string: the method to be used, an empty string
          indicates '"internal"'. }
  \item{tempoutdir}{ is the name of the directory in which the files downloaded on the ArrayExpress repository will be extracted.}
}
\value{
  \code{ Success is indicated by returning the directory in which the
    files have been extracted. If it fails, it returns an empty
    character string.}
}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\keyword{datasets}
