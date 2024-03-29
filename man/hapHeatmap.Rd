% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphic_functions.R
\name{hapHeatmap}
\alias{hapHeatmap}
\title{Graphical output of alleles division by chromosome}
\usage{
hapHeatmap(
  hap_table,
  chain = c("IGH", "IGK", "IGL", "TRB", "TRA"),
  genes_order = NULL,
  removeIGH = TRUE,
  lk_cutoff = 1,
  mark_low_lk = TRUE,
  size_annot = 1.5,
  color_y = NULL,
  order_subject = NULL,
  file = NULL,
  size_text = NULL,
  ylabel_size = 1
)
}
\arguments{
\item{hap_table}{haplotype summary table. See details.}

\item{chain}{the IG chain: IGH,IGK,IGL. Default is IGH.}

\item{genes_order}{A vector of the genes by the desired order. Default is by GENE.loc}

\item{removeIGH}{if TRUE, 'IGH'\'IGK'\'IGL'\'TRB' prefix is removed from gene names.}

\item{lk_cutoff}{the lK cutoff value to be considered low for texture layer. Default is lK<1.}

\item{mark_low_lk}{if TRUE, a texture is add for low lK values. Default is TRUE.}

\item{size_annot}{size of bottom annotation text. Default is 1.5 .}

\item{color_y}{named list of the colors for y axis labels.}

\item{order_subject}{order subject by a vector.}

\item{file}{file path for rendering the plot to pdf. If non is supplied than the plot is returned as object. Default is NULL.}

\item{size_text}{text size for annotations.}

\item{ylabel_size}{text size for y axis labels.}
}
\value{
A list with the following:

\itemize{
  \item \code{'p'}:        heat-map visualization of the haplotype inference for multiple samples.
  \item \code{'width'}:    Optimal width value for rendering plot.
  \item \code{'height'}:   Optimal width value for rendering plot.
}

When a file is supplied the graph is also rendered to pdf.
}
\description{
The \code{hapHeatmap} function generates a graphical output of the alleles per gene in multiple samples.
}
\details{
A \code{data.frame} created by \code{createFullHaplotype}.
}
\examples{
# Plotting haplotpe heatmap
p <- hapHeatmap(samplesHaplotype)
p$p
}
