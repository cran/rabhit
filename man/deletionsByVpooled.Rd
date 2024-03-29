% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{deletionsByVpooled}
\alias{deletionsByVpooled}
\title{Single chromosomal D or J gene deletions inferred by the V pooled method}
\usage{
deletionsByVpooled(
  clip_db,
  chain = c("IGH", "IGK", "IGL"),
  deletion_col = c("d_call", "j_call"),
  count_thresh = 50,
  deleted_genes = "",
  min_minor_fraction = 0.3,
  kThreshDel = 3,
  nonReliable_Vgenes = c()
)
}
\arguments{
\item{clip_db}{a \code{data.frame} in AIRR format. See details.}

\item{chain}{the IG chain: IGH,IGK,IGL. Default is IGH.}

\item{deletion_col}{a vector of column names for which single chromosome deletions should be inferred. Default is j_call and d_call.}

\item{count_thresh}{integer, the minimun number of sequences mapped to a specific V gene to be included in the V pooled inference.}

\item{deleted_genes}{double chromosome deletion summary table. A \code{data.frame} created by \code{deletionsByBinom}.}

\item{min_minor_fraction}{the minimum minor allele fraction to be used as an anchor gene. Default is 0.3}

\item{kThreshDel}{the minimum lK (log10 of the Bayes factor) to call a deletion. Default is 3.}

\item{nonReliable_Vgenes}{a list of known non reliable gene assignments. A \code{list} created by \code{nonReliableVGenes}.}
}
\value{
A \code{data.frame}, in which each row is the single chomosome deletion inference of a gene.

The output containes the following columns:
\itemize{
 \item \code{subject}:       the subject name.
 \item \code{gene}:          the gene call
 \item \code{deletion}:      chromosome deletions inferred. Encoded 1 for deletion and 0 for no deletion.
 \item \code{k}:             the Bayesian factor value for the deletion inference.
 \item \code{counts}:        the appereance count of the gene on each chromosome, the counts are seperated by a comma.
}
}
\description{
The \code{deletionsByVpooled} function inferes single chromosomal deletion for D and J gene .
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
\examples{
\donttest{
data(samples_db)

# Infering V pooled deletions
del_db <- deletionsByVpooled(samples_db)
head(del_db)
}
}
