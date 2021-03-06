% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{read_MAF}
\alias{read_MAF}
\title{Create mutational matrix from MAF file}
\usage{
read_MAF(path, ...)
}
\arguments{
\item{path}{path to MAF file}

\item{...}{other maftools::mutCountMatrix arguments}
}
\value{
the mutational (sparse) matrix associated to the MAF file
}
\description{
Read a MAF (Mutation Annotation Format) file and create
a Mutational Matrix combining gene and sample IDs.
}
\examples{
read_MAF(system.file("extdata", "paac_jhu_2014_500.maf", package = "CIMICE", mustWork = TRUE))

}
