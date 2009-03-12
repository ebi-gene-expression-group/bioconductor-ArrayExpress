\name{extract.zip}
\alias{extract.zip}
\docType{data}
\title{ Unzip archives in a specified directory }
\description{
  \code{extract.zip} extracts the files from a .zip archive in a
  specific directory.
}
\usage{
extract.zip(file, unzip = getOption("unzip"), extract_path = dirname(file))
}
\arguments{
  \item{file}{ A file name. }
  \item{unzip}{ character string: the method to be used, an empty string
          indicates 'internal'. }
  \item{extract_path}{ a path to define where the files are to be extracted. }
}
\value{
  \code{ Success is indicated by returning the directory in which the
    files have been extracted. If it fails, it returns an empty
    character string. }
}
\author{
  Audrey Kauffmann 
  Maintainer: <audrey@ebi.ac.uk>  
}
\keyword{datasets}
