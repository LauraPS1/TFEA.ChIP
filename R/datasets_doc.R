#' Metadata data frame
#'
#' Used to run examples. Data frame containing metadata information
#' for the ChIP-Seq GSM2390643.
#' Fields in the data frame:
#' \itemize{
#'   \item Name: Name of the file.
#'   \item Accession: Accession ID of the experiment.
#'   \item Cell: Cell line or tissue.
#'   \item 'Cell Type': More information about the cells.
#'   \item Treatment
#'   \item Antibody
#'   \item TF: Transcription factor tested in the ChIP-Seq experiment.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ARNT.metadata
#' @usage data("ARNT.metadata")
#' @format a data frame of one row and 7 variables.
"ARNT.metadata"

#' ChIP-Seq dataset
#'
#' Used to run examples. Data frame containing peak information
#' from the ChIP-Seq GSM2390643.
#' Fields in the data frame:
#' \itemize{
#'   \item Name: Name of the file.
#'   \item chr: Chromosome, factor
#'   \item start: Start coordinate for each peak
#'   \item end: End coordinate for each peak
#'   \item X.10.log.pvalue.: log10(p-Value) for each peak.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ARNT.peaks.bed
#' @usage data("ARNT.peaks.bed")
#' @format a data frame of 2140 rows and 4 variables.
"ARNT.peaks.bed"

#' DHS databse
#'
#' Used to run examples. Part of a DHS database storing 76 sites
#' for the human genome in GenomicRanges format.
#'
#' @docType data
#' @keywords datasets
#' @name DnaseHS_db
#' @usage data("DnaseHS_db")
#' @format GenomicRanges object with 76 elements
"DnaseHS_db"

#' List of Entrez Gene IDs
#'
#' Used to run examples. Array of 2754 Entrez Gene IDs extracted
#' from an RNA-Seq experiment sorted by log(Fold Change).
#'
#' @docType data
#' @keywords datasets
#' @name Entrez.gene.IDs
#' @usage data("Entrez.gene.IDs")
#' @format Array of 2754 Entrez Gene IDs.
"Entrez.gene.IDs"

#' List of Entrez Gene IDs
#'
#' Used to run examples. Array of 342 Entrez Gene IDs extracted from
#' upregulated genes in an RNA-Seq experiment.
#'
#' @docType data
#' @keywords datasets
#' @name Genes.Upreg
#' @usage data("Genes.Upreg")
#' @format Array of 2754 Entrez Gene IDs.
"Genes.Upreg"

#' List of one ChIP-Seq dataset
#'
#' Used to run examples. List of part of one ChIP-Seq dataset (from
#' wgEncodeEH002402) in GenomicRanges format with 50 peaks.
#'
#' @docType data
#' @keywords datasets
#' @name gr.list
#' @usage data("gr.list")
#' @format List of one ChIP-Seq dataset.
"gr.list"

#' Output of the function GSEA.run from the TFEA.ChIP package
#'
#' Used to run examples. Output of the function GSEA.run from the
#' TFEA.ChIP package, contains an enrichment table and two lists,
#' one storing runnign enrichment scores and the other,
#' matches/missmatches along a gene list.
#'
#' @docType data
#' @keywords datasets
#' @name GSEA.result
#' @usage data("GSEA.result")
#' @format list of three elements, an erihcment table (data frame), and two list of arrays.
"GSEA.result"

#' RNA-Seq experiment
#'
#' A data frame containing information of of an RNA-Seq experiment
#' on newly transcripted RNA in HUVEC cells during two conditions,
#' 8h of normoxia and 8h of hypoxia (deposited at GEO as GSE89831).
#' The data frame contains the following fields:
#' \itemize{
#'   \item Gene: Gene Symbol for each gene analyzed.
#'   \item Log2FoldChange: base 2 logarithm of the fold change on RNA transcription for a given gene between the two conditions.
#'   \item pvalue
#'   \item padj: p-value adjusted via FDR.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name hypoxia
#' @usage data("hypoxia")
#' @format a data frame of 17527 observations of 4 variables.
"hypoxia"

#' RNA-Seq experiment
#'
#' A DESeqResults objetc containing information of of an RNA-Seq
#' experiment on newly transcripted RNA in HUVEC cells during
#' two conditions, 8h of normoxia and 8h of hypoxia (deposited
#' at GEO as GSE89831).
#'
#' @docType data
#' @keywords datasets
#' @name hypoxia_DESeq
#' @format a DESeqResults objtec
"hypoxia_DESeq"

#' List of Entrez Gene IDs
#'
#' Used to run examples. Array of 2754 log2(Fold Change) values
#' extracted from an RNA-Seq experiment.
#'
#' @docType data
#' @keywords datasets
#' @name log2.FC
#' @usage data("log2.FC")
#' @format Array of 2754 log2(Fold Change) values.
"log2.FC"

#' TF-Gene List
#'
#' This dataset contains two elements: the first is "Gene Keys," which
#' includes all the Entrez IDs, and the second is "ChIP Targets," a list
#' containing information from multiple ChIP-Seq experiments. Each entry
#' in this list contains the indices of the Entrez IDs from the first
#' element that are associated with the peaks of that specific ChIP-Seq.
#'
#' @docType data
#' @keywords datasets
#' @name ChIPDB
#' @usage data("ChIPDB")
#' @format A list with two elements:
#' Gene Keys, a vector of Entrez IDs
#' ChIP Targets, a list containing indices of Entrez IDs associated with peaks from various ChIP-Seq experiments.
#'
"ChIPDB"


#' TF-gene binding DB metadata
#'
#' A data frame containing information about the ChIP-Seq experiments
#' used to build the TF-gene binding DB.
#' Fields in the data frame:
#' \itemize{
#'   \item Accession: Accession ID of the experiment.
#'   \item Cell: Cell line or tissue.
#'   \item 'Cell Type': More information about the cells.
#'   \item Treatment
#'   \item Antibody
#'   \item TF: Transcription factor tested in the ChIP-Seq experiment.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name MetaData
#' @usage data("MetaData")
#' @format A data frame of 1060 observations of 6 variables
"MetaData"


