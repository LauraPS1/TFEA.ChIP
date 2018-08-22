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

#' List of contingency matrix
#'
#' Used to run examples. List of 10 contingency matrix, output of
#' the function "contingency_matrix" from the TFEA.ChIP package.
#'
#' @docType data
#' @keywords datasets
#' @name CM_list
#' @usage data("CM_list")
#' @format a list of 10 contingency matrix.
"CM_list"

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

#' TF-gene binding binary matrix
#'
#' Its rows correspond to all the human genes in the Known Gene database,
#' and its columns, to every ChIP-Seq experiment in the database. The values
#' are 1 – if the ChIP-Seq has a peak assigned to that gene – or 0 –
#' if it hasn’t –.
#'
#' @docType data
#' @keywords datasets
#' @name Mat01
#' @usage data("Mat01")
#' @format a matrix of 1154 columns and 23056 rows
"Mat01"

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
#' @format A data frame of 1154 observations of 6 variables
"MetaData"

#' Data frame, output from the function getCMstats from the TFEA.ChIP package
#'
#' Used to run examples. Output of the function getCMstats from the
#' TFEA.ChIP package, is a data frame storing the following fields:
#' \itemize{
#'   \item Accession: GEO or Encode accession ID for each ChIP-Seq dataset.
#'   \item Cell: cell type on which the ChIP-Seq experiment was performed
#'   \item Treatment: treatment used on the cells
#'   \item TF: Transcription Factor tested.
#'   \item p.value: raw p-value of the Fisher test performed on a contingency matrix for each ChIP-Seq experiment.
#'   \item OR: Odds Ratio on the contingency matrix done for each ChIP-Seq experiment.
#'   \item log2.OR
#'   \item adj.p.value: p-value adjusted by FDR
#'   \item log10.adj.pVal
#'   \item distance: euclidean distance from (log10.adj.pval, log2.OR) to the coordinates origin.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name stat_mat
#' @usage data("stat_mat")
#' @format a data frame of 10 rows and 6 variables
"stat_mat"

#' TFBS database for 3 ChIP-Seq datasets.
#'
#' Used to run examples. Output of the function GR2tfbs_db from the
#' TFEA.ChIP package. Contains a list of three vectors of Entrez Gene
#' IDs assoctiated to three ChIP-Seq experiments
#'
#' @docType data
#' @keywords datasets
#' @name tfbs.database
#' @usage data("tfbs.database")
#' @format a data frame of 10 rows and 6 variables
"tfbs.database"

