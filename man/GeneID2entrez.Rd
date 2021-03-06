% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/All-functions.R
\name{GeneID2entrez}
\alias{GeneID2entrez}
\title{Translates gene IDs from Gene Symbol or Ensemble ID to Entrez ID.}
\usage{
GeneID2entrez(gene.IDs, return.Matrix = FALSE, mode = "h2h")
}
\arguments{
\item{gene.IDs}{Array of Gene Symbols or Ensemble Gene IDs.}

\item{return.Matrix}{T/F. When TRUE, the function returns a matrix[n,2],
one column with the gene symbols or Ensemble IDs, another with their
respective Entrez IDs.}

\item{mode}{Specify the organism used: 'h2h' for homo sapiens gene IDs,
'm2m' for mouse gene IDs, or 'm2h' to get the corresponding human gene
IDs from a mouse input.}
}
\value{
Vector or matrix containing the Entrez IDs(or NA) corresponding
to every element of the input.
}
\description{
Translates mouse or human gene IDs from Gene Symbol or
Ensemble Gene ID to Entrez Gene ID using the IDs approved by HGNC.
When translating from Gene Symbol, keep in mind that many genes have
been given more than one symbol through the years. This function will
return the Entrez ID corresponding to the currently approved symbols
if they exist, otherwise NA is returned. In addition some genes might
map to more than one Entrez ID, in this case gene is assigned to the
first match and a warning is displayed.
}
\examples{
GeneID2entrez(c('TNMD','DPM1','SCYL3','FGR','CFH','FUCA2','GCLC'))
}
