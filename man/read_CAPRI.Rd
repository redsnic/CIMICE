% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{read_CAPRI}
\alias{read_CAPRI}
\title{Read a "CAPRI" file}
\usage{
read_CAPRI(filepath)
}
\arguments{
\item{filepath}{path to file}
}
\value{
the described mutational matrix as a (sparse) matrix
}
\description{
Read a "CAPRI" formatted file from the file system
}
\examples{
#          "pathToDataset/myDataset.CAPRI"
read_CAPRI(system.file("extdata", "example.CAPRI", package = "CIMICE", mustWork = TRUE))

}
