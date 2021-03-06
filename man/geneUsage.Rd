% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{geneUsage}
\alias{geneUsage}
\title{Double chromosome deletion by relative gene usage}
\usage{
geneUsage(clip_db, chain = c("IGH", "IGK", "IGL"))
}
\arguments{
\item{clip_db}{a \code{data.frame} in AIRR format. See details.}

\item{chain}{the IG chain: IGH,IGK,IGL. Default is IGH.}
}
\value{
A \code{data.frame}, in which each row is the relative gene usage value per individual.

The output containes the following columns:
\itemize{
 \item \code{subject}:       the subject name.
 \item \code{gene}:          the gene call
 \item \code{frac}:          the relative gene usage of the gene
}
}
\description{
The \code{geneUsage} function calculates the relative gene usage.
}
\details{
The function accepts a \code{data.frame} in AIRR format (\url{https://changeo.readthedocs.io/en/stable/standard.html}) containing the following columns:
\itemize{
  \item \code{'subject'}: The subject name
  \item \code{'v_call'}: V allele call(s) (in an IMGT format)
  \item \code{'d_call'}: D allele call(s) (in an IMGT format, only for heavy chains)
  \item \code{'j_call'}: J allele call(s) (in an IMGT format)
}
}
