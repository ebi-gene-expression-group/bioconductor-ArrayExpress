\name{getcolproc}
\alias{getcolproc}
\docType{data}
\title{ Return the possible column names from processed MAGE-TAB files }
\description{
  \code{getcolproc} extracts the column names from processed MAGE-TAB and
  return them. The output is needed to call the function \code{procset}.
}

\usage{
getcolproc(files)
}

\arguments{
  \item{files}{ A list as given from \code{\link[ArrayExpress]{getAE}}
    function. Containing the following elements:  
    \describe{

  \item{procfile}{procfile is the name of the processed MAGE-TAB file to
  be read.} 
  \item{path}{is the name of the directory where to find this file.}
  }
}
}

\seealso{\code{\link[ArrayExpress]{ArrayExpress}},
  \code{\link[ArrayExpress]{queryAE}},
  \code{\link[ArrayExpress]{getAE}},
  \code{\link[ArrayExpress]{procset}}} 
\author{
  Audrey Kauffmann
  
  Maintainer: <audrey@ebi.ac.uk>  
}
%\examples{
%# An example can be found in the help of the getAE function.
%}
\keyword{datasets}
