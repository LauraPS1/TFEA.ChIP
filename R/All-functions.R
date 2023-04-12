################### IMPORTS ###################
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges GRanges distanceToNearest mcols
#' @importFrom IRanges IRanges
#' @importFrom biomaRt select useMart getLDS getBM
#' @importFrom dplyr '%>%' arrange
#' @importFrom grDevices colorRamp
#' @importFrom stats fisher.test p.adjust
#' @importFrom utils data setTxtProgressBar txtProgressBar
#' @importFrom R.utils withTimeout
#' @importFrom methods is
#' @import org.Hs.eg.db
################### FUNCTIONS #################

##### Create a ChIP-Gene data base ####
txt2GR <- function( fileTable, format, fileMetaData, alpha = NULL ) {

    #' @title Function to filter a ChIP-Seq input.
    #' @description Function to filter a ChIP-Seq output (in .narrowpeak or
    #' MACS's peaks.bed formats) and then store the peak coordinates in a
    #' GenomicRanges object, associated to its metadata.
    #' @param fileTable data frame from a txt/tsv/bed file
    #' @param format 'narrowpeak', 'macs1.4' or 'macs2'.
    #' narrowPeak fields:
    #' 'chrom','chromStart','chromEnd','name','score','strand','signalValue',
    #' 'pValue','qValue','peak'
    #' macs1.4 fields:
    #' 'chrom','chromStart','chromEnd','name','-10*log10(p-value)'
    #' macs2 fields:
    #' 'chrom','chromStart','chromEnd','name','-log10(p-value)'
    #' @param fileMetaData Data frame/matrix/array contaning the following
    #' fields: 'Name','Accession','Cell','Cell Type','Treatment','Antibody',
    #' 'TF'.
    #' @param alpha max p-value to consider ChIPseq peaks as significant and
    #' include them in the database. By default alpha is 0.05 for narrow peak
    #' files and 1e-05 for MACS files
    #' @return The function returns a GR object generated from the ChIP-Seq
    #' dataset input.
    #' @export txt2GR
    #' @examples
    #' data('ARNT.peaks.bed','ARNT.metadata',package = 'TFEA.ChIP')
    #' ARNT.gr<-txt2GR(ARNT.peaks.bed,'macs1.4',ARNT.metadata)

    stopifnot(format %in% c("narrowpeak","narrowPeak","macs1.4",
                            "macs2", "MACS1.4", "MACS2"))
    stopifnot((is.data.frame(fileMetaData) || is.matrix(fileMetaData) ||
        is.array(fileMetaData)))

    # checking fileMetadata has all the fields needed and the correct
    # column names and column order.
    columnNames <- c("Name", "Accession", "Cell", "Cell.Type",
        "Treatment", "Antibody", "TF")


    fileMetaData <- as.data.frame( fileMetaData )

    if ( dim( fileMetaData )[2] == 7) {

        if (FALSE %in% (columnNames %in% colnames(fileMetaData))) {
            stop("fileMetaData format error: 'fileMetaData' must be a",
                " data frame/matrix/array with 7 atributes: 'Name',",
                "'Accession', 'Cell', 'Cell.Type','Treatment',",
                "'Antibody','TF'")
        } else {
            fileMetaData <- fileMetaData[columnNames]
        }
    } else {
        stop("fileMetaData format error: 'fileMetaData' must be a",
            " data frame/matrix/array with 7 atributes: 'Name',",
            "'Accession', 'Cell', 'Cell.Type','Treatment',",
            "'Antibody','TF'")
    }

    format <- tolower(format)

    if (format == "narrowpeak") {

        if (fileTable[1, 8] == -1 & fileTable[1, 9] == -1) {
            # If there's no p-value column
            warning("The ChIP-Seq input file ", fileMetaData$Name,
                " does not include p-value or Q-value for each peak.",
                " Please, make sure the peaks in the input file have been",
                " previously filtered according to their significance")
            fileTable <- fileTable[, 1:3]
            Stat <- "no score"
            fileTable$score = rep(NA, dim(fileTable)[1])
            colnames(fileTable)[1:3] <- c("chr", "start", "end")

        } else if (fileTable[1, 8] == -1) {
            # If the score table has a q-value column
            fileTable <- fileTable[, c(1, 2, 3, 9)]
            colnames(fileTable) <- c("chr", "start", "end", "score")
            Stat <- "log10(p-Value)"

            if (is.null(alpha)) {
                valLimit <- 1.3
            } else {
                valLimit <- (-log10(alpha))
            }
            fileTable <- fileTable[fileTable$score > valLimit, ]

        } else if (fileTable[1, 9] == -1) {
            # If the score table has a -log10(p-value) column
            fileTable <- fileTable[, c(1, 2, 3, 8)]
            colnames(fileTable) <- c("chr", "start", "end", "score")
            fileTable$score <- 10 ^ (-1 * (fileTable$score ) )
            fileTable$score <- p.adjust( fileTable$score, "fdr" ) # adjust p-values using Benjamini & Hochberg correction.
            Stat <- "corrected p-Value"
            if (is.null(alpha)) {
                valLimit <- 0.05
            } else {
                valLimit <- alpha
            }
            fileTable <- fileTable[fileTable$score < valLimit, ]

        } else {
            # If the dataset has both p-value and 10*-log(qval) columns.
            fileTable <- fileTable[, c(1, 2, 3, 9)]
            colnames(fileTable) <- c("chr", "start", "end", "score")
            Stat <- "corrected p-Value"

            if (is.null(alpha)) {
                valLimit <- 1.3
            } else {
                valLimit <- (-log10(alpha))
            }
            fileTable <- fileTable[fileTable$score > valLimit, ]
        }
        if( dim( fileTable)[1] > 0){
            fileMetaData <- c(fileMetaData, Stat)
            MDframe <- as.data.frame(lapply(fileMetaData, rep, dim(fileTable)[1]))
            colnames(MDframe) <- c("Name", "Accession", "Cell", "Cell Type",
                "Treatment", "Antibody", "TF", "Score Type")

            gr <- GenomicRanges::GRanges(seqnames = fileTable$chr,
                ranges = IRanges::IRanges(fileTable$start, end = fileTable$end),
                score = fileTable$score, mcols = MDframe)
            return(gr)
        } else {
            warning( "File ", fileMetaData$Accession,
                     " has no significant peaks after correcting p-values.")
            return(NULL)
        }

    } else if (format %in% c("macs1.4", "macs2") ) {

        if (is.null(alpha)) {
            valLimit <- 0.05
        } else {
            valLimit <- alpha
        }

        if (dim(fileTable)[2] == 5) {
            fileTable <- fileTable[, c(1, 2, 3, 5)]
            colnames(fileTable) <- c("chr", "start", "end", "score")

            if (format == "macs1.4") {
                fileTable$score <- 10 ^ (-1 * (fileTable$score / 10 ) ) # -10 * log10(p-value)
            } else if (format == "macs2"){
                fileTable$score <- 10 ^ (-1 * (fileTable$score ) ) # -log10(p-value)
            }
            fileTable$score <- p.adjust( fileTable$score, "fdr" ) # adjust p-values using Benjamini & Hochberg correction.

            fileTable <- fileTable[fileTable$score < valLimit, ]
            Stat <- "corrected p-Value"

        } else if (dim(fileTable)[2] == 4 & is.character(fileTable[1, 4])) {
            # if the 4th column consists of peak names
            warning("The ChIP-Seq input file does not include p-value or",
                " Q-value for each peak. Please, make sure the peaks",
                " in the input file have been previously filtered",
                " according to their significance")
            fileTable <- fileTable[, 1:3]
            Stat <- "no score"
            fileTable$score = rep(NA, dim(fileTable)[1])
            colnames(fileTable)[1:3] <- c("chr", "start", "end")

        } else if (dim(fileTable)[2] == 4 & !is.character(fileTable[1,
            4])) {
            # if the 4th column consists of p-values
            colnames(fileTable) <- c("chr", "start", "end", "score")

            if (format == "macs1.4") {
                fileTable$score <- 10 ^ (-1 * (fileTable$score / 10 ) ) # -10 * log10(p-value)
            } else if (format == "macs2"){
                fileTable$score <- 10 ^ (-1 * (fileTable$score ) ) # -log10(p-value)
            }

            fileTable$score <- p.adjust( fileTable$score, "fdr" ) # adjust p-values using Benjamini & Hochberg correction.
            fileTable <- fileTable[fileTable$score < valLimit, ]
            Stat <- "corrected p-Value"
        }

        if( dim( fileTable)[1] > 0){

            fileMetaData <- c(fileMetaData, Stat)
            MDframe <- as.data.frame(lapply(fileMetaData, rep, dim(fileTable)[1]))
            colnames(MDframe) <- c("Name", "Accession", "Cell", "Cell Type",
                "Treatment", "Antibody", "TF", "Score Type")

            gr <- GenomicRanges::GRanges(seqnames = fileTable$chr,
                ranges = IRanges::IRanges(fileTable$start, end = fileTable$end),
                score = fileTable$score, mcols = MDframe)
            return(gr)
        } else {
            warning( "File ", fileMetaData$Accession,
                     " has no significant peaks after correcting p-values.")
            return( NULL )
        }

    } else {
        stop("format error: variable 'format' must be",
             " either 'narrowpeak' or 'macs'. ")
    }
}

makeChIPGeneDB <- function( Ref.db, gr.list, distanceMargin = 10,
                            min.Targets = 10 ){
    #' @title Make a ChIP - target database
    #' @description makeChIPGeneDB generates a ChIP-seq - target database
    #' through the association of ChIP-Seq peak coordinates (provided
    #' as a GenomicRange object) to overlapping genes or gene-associated
    #' genomic regions (Ref.db).
    #' @param Ref.db GenomicRanges object containing a database of reference
    #' elements (either Genes or gene-associated regions) including a gene_id
    #' metacolumn
    #' @param gr.list List of GR objects containing ChIP-seq peak coordinates
    #' (output of txt2GR).
    #' @param distanceMargin Maximum distance allowed between a gene or 
    #' regulatory element to assign a gene to a ChIP-seq peak. Set to 10 bases
    #' by default.
    #' @param min.Targets Minimum number of putative targets per ChIP-seq in
    #' gr.list. ChIPs with fewer targets will be discarded.
    #' regulatory element to assign a gene to a ChIP-seq peak. Set to 10 bases
    #' by default.
    #' @return List containing two elements:
    #'    - Gene Keys: vector of gene IDs
    #'    - ChIP Targets: list of vectors, one per element in gr.list, 
    #'      containing the putative targets assigned. Each target is coded as
    #'      its position in the vector 'Gene Keys'.
    #' @export makeChIPGeneDB
    #' @examples
    #' data( 'DnaseHS_db','gr.list', package = 'TFEA.ChIP' )
    #' makeChIPGeneDB( DnaseHS_db, gr.list )
  
  
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("S4Vectors package needed for this function to work. ",
         "Please install it.", call. = FALSE)
  }
  
  referenceIDs <- sort(unique( Ref.db$gene_id ))
  
  ChIPtarget_list <- lapply(
    gr.list,
    function( gr.list, dm, Ref.db, IDs, min.Targets ,i ){
      
      nearest_index <- suppressWarnings(
        GenomicRanges::findOverlaps( i, Ref.db, maxgap = dm ))
      inSubject <- S4Vectors::subjectHits( nearest_index )
      
      assigned_genes <- unique( Ref.db[ inSubject ]$gene_id )
      
      # in case any ChIP-Seq dataset does not have enough targets
      if ( length( assigned_genes ) < min.Targets ) {
        NULL
      } else {
        sort( match( assigned_genes, IDs ) )
      }
    },
    gr.list = gr.list,
    dm = distanceMargin,
    Ref.db = Ref.db, 
    IDs = referenceIDs,
    min.Targets = min.Targets
  )
  
  # naming elements
  if( !is.null( names( gr.list ))){
    names( ChIPtarget_list ) <- names( gr.list )
  }else {
    # if gr.list is not named, search for an 'Accession' in GR metadata
    # or create names.
    ChIP_acc_col <- grep(
      "accession",
      tolower( colnames( mcols( gr.list[[ 1 ]] ) ) ) )[1]
    
    if ( length( ChIP_acc_col > 0 ) ){
      ChIP_acc_col <- ChIP_acc_col[ 1 ]
      cat("Using '", ChIP_acc_col, "' as ChIP-seq identifiers." )
      
      chipNames <- sapply(
        gr.list,
        function(i){
          mcols( i )@listData[[ ChIP_acc_col ]][ 1 ]
        }
      )
      names( ChIPtarget_list ) <- chipNames 
      
    } else {
      names( ChIPtarget_list ) <- paste0("ChIPseq_", 1:length(gr.list ) )
    }
  }
  
  ChIPtarget_list[ sapply( ChIPtarget_list, is.null ) ] <- NULL
  
  return( list(
    "Gene Keys" = referenceIDs,
    "ChIP Targets" = ChIPtarget_list
  ) )
  
}

matrixDB_to_listDB <- function( Mat01 ){
    #' @title Re-formatting ChIP-Gene database
    #' @description Function to transform a ChIP-gene data base from the former
    #' binary matrix to the current list-based format.
    #' @param Mat01 Matrix[n,m] which rows correspond to all the human
    #' genes that have been assigned an Entrez ID, and its columns, to every
    #' ChIP-Seq experiment in the database. The values are 1 – if the ChIP-Seq
    #' has a peak assigned to that gene – or 0 – if it hasn’t –.
    #' @return List containing two elements:
    #'    - Gene Keys: vector of gene IDs
    #'    - ChIP Targets: list of vectors, one per ChIP-seq experiment in the,
    #'      database, containing the putative targets assigned. Each target is
    #'      coded as its position in the vector 'Gene Keys'.
    #' @export matrixDB_to_listDB
    #' @examples
    #' Mat01 <- matrix( 
    #'     round( runif(9) ), nrow = 3,
    #'     dimnames= list( paste0("Gene ", 1:3), paste0("ChIPseq ", 1:3))  )
    #' matrixDB_to_listDB( Mat01 )
    
    ChIPDB <- list(
      "Gene Keys" = rownames( Mat01 )
    )
    ChIPDB <- c( ChIPDB, "ChIP Targets" = list( sapply(
        seq_len( ncol( Mat01 ) ),
        function(i){
            unname( which( Mat01[, i] == 1 ) )
        }
    )) )
    names( ChIPDB[["ChIP Targets"]] ) <- colnames( Mat01 )
    
    return( ChIPDB )
}

#### Prepare input data and options ####

set_user_data <- function( metadata, ChIPDB ) {
    #' @title Sets the data objects as default.
    #' @description Function to set the data objects provided by the user
    #' as default to the rest of the functions.
    #' @param metadata Data frame/matrix/array contaning the following fields:
    #' 'Name','Accession','Cell','Cell Type','Treatment','Antibody','TF'.
    #' @param ChIPDB List containing two elements:
    #'    - Gene Keys: vector of gene IDs
    #'    - ChIP Targets: list of vectors, one per ChIP-seq experiment in the,
    #'      database, containing the putative targets assigned. Each target is
    #'      coded as its position in the vector 'Gene Keys'.
    #' @return sets the user's metadata table and TFBS matrix as the variables
    #' 'MetaData' and 'ChIPDB', used by the rest of the package.
    #' @export set_user_data
    #' @examples
    #' data( 'MetaData', 'ChIPDB', package='TFEA.ChIP' )
    #' # For this example, we will usethe variables already included in the
    #' # package.
    #' set_user_data( MetaData, ChIPDB )

    pos <- 1
    envir = as.environment(pos)

    assign("MetaData", metadata, envir = envir)
    assign("ChIPDB", ChIPDB, envir = envir)
}

preprocessInputData <- function(inputData, mode = "h2h" ) {
    #' @title Extracts data from a DESeqResults object or a data frame.
    #' @description Function to extract Gene IDs, logFoldChange, and p-val
    #' values from a DESeqResults object or data frame. Gene IDs are
    #' translated to ENTREZ IDs, if possible, and the resultant data frame
    #' is sorted accordint to decreasing log2(Fold Change). Translating
    #' gene IDs from mouse to their equivalent human genes is avaible
    #' using the variable "mode".
    #' @param inputData DESeqResults object or data frame. In all cases
    #' must include gene IDs. Data frame inputs should include 'pvalue' and
    #' 'log2FoldChange' as well.
    #' @param  mode Specify the organism used: 'h2h' for homo sapiens gene IDs,
    #' 'm2m' for mouse gene IDs, or 'm2h' to get the corresponding human gene
    #' IDs from a mouse input.
    #' @return A table containing Entrez Gene IDs, LogFoldChange and p-val
    #' values (both raw p-value and fdr adjusted p-value), sorted by
    #' log2FoldChange.
    #' @export preprocessInputData
    #' @examples
    #' data('hypoxia_DESeq',package='TFEA.ChIP')
    #' preprocessInputData( hypoxia_DESeq )

    if ( methods::is(inputData, "DESeqResults") ){
        # Extracting data from a DESeqResults object
        if (!requireNamespace("DESeq2", quietly = TRUE)) {
            stop("DESeq2 package needed for this function to work. ",
                "Please install it.", call. = FALSE)
        }
        g <- inputData@rownames
        inputData <- as.data.frame( inputData@listData )
        rownames( inputData ) <- g
        rm(g)

        # check the gene ids and translate if needed
        if ( ! all( grepl("^\\d*$", rownames(inputData) ) ) | mode == "m2h" ) {
            genes <- suppressMessages( GeneID2entrez(
                gene.IDs = rownames( inputData ),
                mode,
                return.Matrix = TRUE))

            if ( mode %in% c("h2h","m2m") ){
                genes <- genes[!is.na(genes$ENTREZ.ID), ]
                inputData <- inputData[rownames(inputData) %in% genes$GENE.ID, ]
                genes <- genes$ENTREZ.ID
            } else {
                genes <- genes[!is.na(genes$human.gene.ID), ]
                inputData <- inputData[rownames(inputData) %in% genes$mouse.gene.ID, ]
                genes <- genes$human.gene.ID
            }

        } else {
            genes <- rownames(inputData)
        }
        # get the rest of the variables
        log2FoldChange <- inputData[["log2FoldChange"]]
        pvalue <- inputData[["pvalue"]]
        pval.adj <- inputData[["padj"]]

        Table <- data.frame(Genes = genes, log2FoldChange = log2FoldChange,
            pvalue = pvalue, pval.adj = pval.adj)
        Table$Genes <- as.character(Table$Genes)
        Table <- Table[!is.na(Table$log2FoldChange), ]
        Table <- Table[order(Table$log2FoldChange, decreasing = TRUE), ]
        rownames(Table) <- NULL
        return(Table)

    } else if ( methods::is(inputData, "data.frame") ) {
        # Extracting data from a data frame.
        # Checkig if all the necessary columns are present
        if (FALSE %in% (c("Genes", "pvalue", "log2FoldChange") %in%
            colnames(inputData)) & FALSE %in% (c("Genes", "pval.adj",
            "log2FoldChange") %in% colnames(inputData))) {
            stop("The necessary atributes can't be found in input data frame",
                ". Input data must include: 'Genes', 'log2FoldChange',and ",
                "'pvalue' or 'pval.adj'", call. = FALSE)
        }
        # If there's not an adjusted p-value column
        if (!("pval.adj" %in% colnames(inputData))) {
            inputData$pval.adj <- p.adjust(inputData$pvalue,
                "fdr")
        }
        # If Gene IDs aren't in Entrez Gene ID format or come from mouse genes.
        if ( ! all( grepl("^\\d*$", inputData$Genes)) | mode == "m2h" ) {
            if( mode== "h2h" ){ inputData$Genes <- toupper( inputData$Genes ) }
            inputData$Genes <- trimws( inputData$Genes )
            genes <- suppressMessages( GeneID2entrez(
                gene.IDs = inputData$Genes,
                mode,
                return.Matrix = TRUE))

            if (  mode %in% c("h2h","m2m")  ){
                genes <- genes[!is.na(genes$ENTREZ.ID), ]
                inputData <- inputData[ inputData$Genes %in% genes$GENE.ID, ]
                inputData$Genes <- genes$ENTREZ.ID
            }else{
                genes <- genes[!is.na(genes$human.gene.ID), ]
                inputData <- inputData[ inputData$Genes %in% genes$mouse.gene.ID, ]
                inputData$Genes <- genes$human.gene.ID
            }
        }
        # sorting according to log2(FoldChange)
        inputData <- inputData[order(inputData$log2FoldChange,
            decreasing = TRUE), ]
        return(inputData)

    } else {
        stop("preprocessInputData requires a DESeqResults object or ",
            "a data frame as input.", call. = FALSE)}
}

Select_genes <- function(GeneExpression_df, max_pval = 0.05,
                         min_pval = 0, max_LFC = Inf, min_LFC = -Inf) {
    #' @title Extracts genes according to logFoldChange and p-val limits
    #' @description Function to extract Gene IDs from a dataframe according
    #' to the established limits for log2(FoldChange) and p-value.
    #' If possible, the function will use the adjusted p-value column.
    #' @param GeneExpression_df A data frame with the folowing fields:
    #' 'Gene', 'pvalue' or 'pval.adj', 'log2FoldChange'.
    #' @param max_pval maximum p-value allowed, 0.05 by default.
    #' @param min_pval minimum p-value allowed, 0 by default.
    #' @param max_LFC maximum log2(FoldChange) allowed.
    #' @param min_LFC minimum log2(FoldChange) allowed.
    #' @return A vector of gene IDs.
    #' @export Select_genes
    #' @examples
    #' data('hypoxia',package='TFEA.ChIP')
    #' Select_genes(hypoxia)
  
  # Checking input variables
  if (max_pval < min_pval) {
    stop("'max_pval' has to be greater than 'min_pval'. ", call. = FALSE)
  }
  if (max_LFC < min_LFC) {
    stop("'max_LFC' has to be greater than 'min_LFC'.", call. = FALSE)
  }
  if( ! "log2FoldChange" %in% colnames( GeneExpression_df ) ) {
    stop("No 'log2FoldChange' field in the input data. ")
  } 
  if( ! any(c( "pvalue","pval.adj" ) %in%  colnames( GeneExpression_df ) ) ){
    stop("No 'pvalue' or 'pval.adj' field in the input data. ")
  }
  
  # Selecting by p-value
  if ("pval.adj" %in% colnames( GeneExpression_df ) ) {
    g_pv <- GeneExpression_df$pval.adj >= min_pval &
      GeneExpression_df$pval.adj <= max_pval
  } else if ("pvalue" %in% colnames(GeneExpression_df)) {
    g_pv <- GeneExpression_df$pvalue >= min_pval &
      GeneExpression_df$pvalue <= max_pval
  }
  
  # Selecting by log2(FoldChange) 
  g_lfc <- GeneExpression_df$log2FoldChange >= min_LFC &
    GeneExpression_df$log2FoldChange <= max_LFC
  
  GeneExpression_df$Gene[ g_pv & g_lfc ]
}

GeneID2entrez <- function(gene.IDs, return.Matrix = FALSE, mode = "h2h") {
    
    #' @title Translates gene IDs from Gene Symbol or Ensemble ID to Entrez ID.
    #' @description Translates mouse or human gene IDs from Gene Symbol or
    #' Ensemble Gene ID to Entrez Gene ID using the IDs approved by HGNC.
    #' When translating from Gene Symbol, keep in mind that many genes have
    #' been given more than one symbol through the years. This function will
    #' return the Entrez ID corresponding to the currently approved symbols
    #' if they exist, otherwise NA is returned. In addition some genes might
    #' map to more than one Entrez ID, in this case gene is assigned to the
    #' first match and a warning is displayed.
    #' @param gene.IDs Array of Gene Symbols or Ensemble Gene IDs.
    #' @param return.Matrix T/F. When TRUE, the function returns a matrix[n,2],
    #' one column with the gene symbols or Ensemble IDs, another with their
    #' respective Entrez IDs.
    #' @param mode Specify the organism used: 'h2h' for homo sapiens gene IDs,
    #' 'm2m' for mouse gene IDs, or 'm2h' to get the corresponding human gene
    #' IDs from a mouse input.
    #' @return Vector or matrix containing the Entrez IDs(or NA) corresponding
    #' to every element of the input.
    #' @export GeneID2entrez
    #' @examples
    #' GeneID2entrez(c('TNMD','DPM1','SCYL3','FGR','CFH','FUCA2','GCLC'))
    
    stopifnot( mode %in% c("h2h", "m2m", "m2h"))
    
    gene.IDs <- gene.IDs[ !is.na( gene.IDs ) ]
    gene.IDs <- trimws( gene.IDs ) # remone any possible white space
    
    if ( mode == 'h2h'){
        
        gene.IDs <- toupper(gene.IDs)  # in case any name is in lowercase.
        # suppressWarnings added to avoid 'select()' returned 1:many
        # mapping between keys and columns
        if ( all( grepl("^ENSG", gene.IDs, perl = TRUE ) ) ) {
            ID.type <- "ENSEMBL"
            suppressMessages(GeneNames <- biomaRt::select(
                org.Hs.eg.db,  gene.IDs,
                c("ENTREZID", "ENSEMBL"), keytype = "ENSEMBL" ))
        } else {
            ID.type <- "SYMBOL"
            suppressMessages(GeneNames <- biomaRt::select(
                org.Hs.eg.db,  gene.IDs,
                c("SYMBOL", "ENTREZID"), keytype = "ALIAS" ))
        }
        
        matched <- match( as.character(gene.IDs), GeneNames[, ID.type] )
        matched_2 <- match( GeneNames[, ID.type], as.character(gene.IDs) )
        
        if ( sum( duplicated( matched_2[ !is.na(matched_2) ] ) ) > 0) {
            warning("Some genes returned 1:many mapping to ENTREZ ID. ",
                    "Genes were assigned the first ENTREZ ID match found.\n",
                    call. = FALSE)
        }
        message("Done! ", sum( !is.na(matched) ), " genes of ",
                length( matched ), " successfully converted.\n")
        
        if ( return.Matrix == TRUE ) {
            if ( sum( is.na(matched) ) > 0 ) {
                message("Couldn't find Entrez IDs for ", sum( is.na(matched) ),
                        " genes (NAs returned instead).\n")
            }
            return( data.frame(
                GENE.ID = gene.IDs,
                ENTREZ.ID = GeneNames[ matched, "ENTREZID"],
                stringsAsFactors = FALSE))
        } else {
            if ( sum( is.na(matched) ) > 0) {
                message("Couldn't find Entrez IDs for ", sum( is.na(matched) ),
                        " genes.\n")
            }
            return( GeneNames[ matched[!is.na(matched)] , "ENTREZID"])
        }
        
    } else if( mode == "m2m" ) {
        
        biomart_test <- tryCatch(
            {R.utils::withTimeout( {tmp <- biomaRt::listMarts()},
                                   timeout = 3, onTimeout = "warning")},
            error = function(w) { return( 0 ) },
            warning = function(w){ return( 0 ) }
        )
        if ( all(biomart_test == 0) ){
            stop( paste0("We are having trouble reaching biomaRt.\n",
                         "Please, try again later."))
        }
        
        mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        GeneNames <- biomaRt::getBM(
            attributes = c("ensembl_gene_id", "mgi_symbol","entrezgene_id"),
            values = "*", mart = mouse)
        
        if ( all( grepl("^ENSM", gene.IDs, perl = TRUE ) ) ) {
            
            ids <- GeneNames[ match( gene.IDs, GeneNames$ensembl_gene_id), c(1,3) ]
            colnames(ids)<- c("GENE.ID", "ENTREZ.ID")
            
        } else {
            
            ids <- GeneNames[ match( gene.IDs, GeneNames$mgi_symbol), c(2,3) ]
            colnames(ids)<- c("GENE.ID", "ENTREZ.ID")
        }
        
        message("Done! ", sum( ! is.na( ids$ENTREZ.ID ) ) ,
                " genes of ", dim( ids )[1], " successfully converted.\n")
        
        if (return.Matrix == TRUE) {
            if ( sum( is.na( ids$ENTREZ.ID )) > 0) {
                message("Couldn't find Entrez IDs for ",
                        sum( is.na( ids$ENTREZ.ID )),
                        " genes (NAs returned instead).\n")
            }
            return( ids )
        } else {
            if ( sum( is.na( ids$ENTREZ.ID )) > 0) {
                message("Couldn't find human Entrez IDs for ",
                        sum( is.na( ids$ENTREZ.ID )), " genes.\n")
            }
            return( ids$ENTREZ.ID[ !is.na( ids$ENTREZ.ID ) ] )
        }
        
        
    } else if( mode == "m2h" ){
        
        biomart_test <- tryCatch(
            {R.utils::withTimeout( {tmp <- biomaRt::listMarts()},
                                   timeout = 3, onTimeout = "warning")},
            error = function(w) { return( 0 ) },
            warning = function(w){ return( 0 ) }
        )
        if( all(biomart_test == 0) ){
            stop( paste0("We are having trouble reaching biomaRt.\n",
                         "Please, try again later."))
        }
        
        
        human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        
        if ( all( grepl("^ENSM", gene.IDs, perl = TRUE ) ) == TRUE) {
            hs_ids = biomaRt::getLDS(
                attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id",
                values = gene.IDs , mart = mouse,
                attributesL = c("entrezgene_id"), martL = human,
                uniqueRows = TRUE )
            hs_ids <- data.frame(
                mouse.gene.ID=gene.IDs,
                human.gene.ID=hs_ids$NCBI.gene.ID[match( gene.IDs, hs_ids$Gene.stable.ID )  ],
                stringsAsFactors = F
            )
            
        } else if ( all( grepl("^\\d*$", gene.IDs, perl = TRUE) ) == TRUE ) {
            hs_ids = biomaRt::getLDS(
                attributes = c("entrezgene_id"), filters = "entrezgene_id",
                values = gene.IDs, mart = mouse,
                attributesL = c("entrezgene_id"), martL = human,
                uniqueRows = TRUE )
            hs_ids <- data.frame(
                mouse.gene.ID=gene.IDs,
                human.gene.ID=hs_ids$NCBI.gene.ID.1[match( gene.IDs, hs_ids$NCBI.gene.ID )  ],
                stringsAsFactors = F
            )
            
        } else {
            hs_ids = biomaRt::getLDS(
                attributes = c("mgi_symbol"), filters = "mgi_symbol",
                values = gene.IDs, mart = mouse,
                attributesL = c("entrezgene_id"), martL = human,
                uniqueRows = TRUE )
            hs_ids <- data.frame(
                mouse.gene.ID=gene.IDs,
                human.gene.ID=hs_ids$NCBI.gene.ID[match( gene.IDs, hs_ids$MGI.symbol )  ],
                stringsAsFactors = F
            )
        }
        
        message("Done! ", sum( ! is.na( hs_ids$human.gene.ID )),
                " genes of ", length( gene.IDs ), " successfully converted.\n")
        
        if (return.Matrix == TRUE) {
            if ( any( is.na( hs_ids$human.gene.ID )) ) {
                message("Couldn't find human Entrez IDs for ",
                        sum(is.na( hs_ids$human.gene.ID )),
                        " genes (NAs returned instead).\n")
            }
            return( hs_ids )
        } else {
            if ( any( is.na( hs_ids$human.gene.ID )) ) {
                message("Couldn't find human Entrez IDs for ",
                        sum(is.na( hs_ids$human.gene.ID )),  " genes.\n")
            }
            return( hs_ids$human.gene.ID[ !is.na( hs_ids$human.gene.ID ) ] )
        }
        
    }
}

get_chip_index <- function(encodeFilter = FALSE, TFfilter = NULL) {

    #' @title Creates df containing accessions of ChIP-Seq datasets and TF.
    #' @description Function to create a data frame containing the ChIP-Seq
    #' dataset accession IDs and the transcription factor tested in each ChIP.
    #' This index is used in functions like “contingency_matrix” and “GSEA_run”
    #' as a filter to select specific ChIPs or transcription factors to run an
    #' analysis.
    #' @param encodeFilter (Optional) If TRUE, only ENCODE ChIP-Seqs are
    #' included in the index.
    #' @param TFfilter (Optional) Transcription factors of interest.
    #' @return Data frame containig the accession ID and TF for every ChIP-Seq
    #' experiment included in the metadata files.
    #' @export get_chip_index
    #' @examples
    #' get_chip_index(encodeFilter = TRUE)
    #' get_chip_index(TFfilter=c('SMAD2','SMAD4'))

    if (!exists("MetaData")) {
        MetaData <- NULL
        data("MetaData", package = "TFEA.ChIP", envir = environment())
    }

    if (is.null(TFfilter)) {
        Index <- dplyr::select(MetaData, "Accession", "TF")
        if (encodeFilter == TRUE) {
            if (any(grepl("^wg.*", Index$Accession))){
                Index <- Index[grepl("^wg", Index$Accession), ]
            } else {
                Index <- Index[grepl("^ENC", Index$Accession), ]
            }
        }
    } else {
        Index <- dplyr::select(MetaData, "Accession", "TF")
        Index <- Index[Index$TF %in% TFfilter, ]
        if (encodeFilter == TRUE) {

            if (any(grepl("^wg.*", Index$Accession))){
                Index <- Index[grepl("^wg", Index$Accession), ]
            } else {
                Index <- Index[grepl("^ENC", Index$Accession), ]
            }
        }
    }

    if ( dim( Index )[1] == 0 ){
        stop("No ChIP-Seq dataset in the database follows ",
            "your conditions.")
    }else{
        return(Index)
    }
}

#### Run TFEA ####

contingency_matrix <- function(test_list, control_list,
                               chip_index = get_chip_index()) {
  
    #' @title Computes 2x2 contingency matrices
    #' @description Function to compute contingency 2x2 matrix by the partition
    #' of the two gene ID lists according to the presence or absence of the
    #' terms in these list in a ChIP-Seq binding database.
    #' @param test_list List of gene Entrez IDs
    #' @param control_list If not provided, all human genes not present in
    #' test_list will be used as control.
    #' @param chip_index Output of the function “get_chip_index”, a data frame
    #' containing accession IDs of ChIPs on the database and the TF each one
    #' tests. If not provided, the whole internal database will be used
    #' @return List of contingency matrices, one CM per element in chip_index
    #' (i.e. per ChIP-seq dataset).
    #' @export contingency_matrix
    #' @examples
    #' data('Genes.Upreg',package = 'TFEA.ChIP')
    #' CM_list_UP <- contingency_matrix(Genes.Upreg)
  
    if ( ! exists( "ChIPDB" )) {
        ChIPDB <- NULL
        data( "ChIPDB", package = "TFEA.ChIP", envir = environment() )
    }
    if( is.matrix( ChIPDB ) ){
        ChIPDB <- matrixDB_to_listDB( ChIPDB )
    }
    
    if ( missing( control_list ) ) {
        # Generating control gene list in case is not provided.
        control_list <- ChIPDB[["Gene Keys"]]  
    }
    control_list <- control_list[ !( control_list %in% test_list ) ]
    
    CM_list <- lapply(
        seq_along(chip_index$Accession),
        
        function( i, g1, g2, accs, ChIPDB ){
            
            acc <- which( names( ChIPDB[[ "ChIP Targets" ]] ) == accs[ i ] )
            chip.targets <- ChIPDB[["Gene Keys"]][
                ChIPDB[[ "ChIP Targets" ]][[ acc ]] ]
            
            pos1 <- sum( g1 %in% chip.targets )
            pos2 <- sum( g2 %in% chip.targets )
            neg1 <- sum( ! g1 %in% chip.targets )
            neg2 <- sum( ! g2 %in% chip.targets )
            
            contMatrix <- cbind(c(pos1, pos2), c(neg1, neg2))
            rownames(contMatrix) <- c("Test", "Control")
            colnames(contMatrix) <- c("Positive", "Negative")
            
            return(contMatrix)
        },
        g1 = test_list,
        g2 = control_list,
        accs = chip_index$Accession,
        ChIPDB = ChIPDB
    )
    
    names(CM_list) <- chip_index$Accession
    
    return(CM_list)
}

getCMstats <- function( CM_list, chip_index = get_chip_index()) {

    #' @title Generate statistical parameters from a contingency_matrix output
    #' @description From a list of contingency matrices, such as the output
    #' from “contingency_matrix”, this function computes a fisher's exact test
    #' for each matrix and generates a data frame that stores accession ID of a
    #' ChIP-Seq experiment, the TF tested in that experiment, the p-value and
    #' the odds ratio resulting from the test.
    #' @param CM_list Output of “contingency_matrix”, a list of
    #' contingency matrices.
    #' @param chip_index Output of the function “get_chip_index”, a data frame
    #' containing accession IDs of ChIPs on the database and the TF each one
    #' tests. If not provided, the whole internal database will be used
    #' @return Data frame containing accession ID of a ChIP-Seq experiment and
    #' its experimental conditions, the TF tested in that experiment, raw and
    #' adjusted p-values, odds-ratio, and euclidean distance.
    #' and FDR-adjusted p-values (-10*log10 adj.pvalue).
    #' @export getCMstats
    #' @examples
    #' data('Genes.Upreg',package = 'TFEA.ChIP')
    #' CM_list_UP <- contingency_matrix( Genes.Upreg )
    #' stats_mat_UP <- getCMstats( CM_list_UP )
    
    
    fishTests <- do.call( rbind, lapply(
        seq_along( CM_list ),
        function(CM_list,i) {
            data.frame( stats::fisher.test(
                CM_list[[ i ]] )[ c("estimate","p.value") ] )
        },
        CM_list = CM_list ) )
    rownames( fishTests ) <- names( CM_list )
    
    chip_index <- chip_index[ match( names(CM_list), chip_index$Accession ),]
    
    statMat <- data.frame(
        Accession = chip_index$Accession, TF = chip_index$TF,
        p.value = fishTests$p.value, OR = fishTests$estimate,
        stringsAsFactors = FALSE)
    
    statMat$log2.OR <- log2(statMat$OR)
    statMat$log2.OR[ abs( statMat$log2.OR ) == Inf ] <- NA

    statMat$adj.p.value <- stats::p.adjust( statMat$p.value, "fdr" )
    statMat$log10.adj.pVal <- ( -1 * ( log10( statMat$adj.p.value ) ) )
    statMat$log10.adj.pVal[ abs( statMat$log2.OR ) == Inf ] <- NA
    
    tmpOR <- statMat$OR
    tmpOR[ statMat$OR == Inf ] <-  max( statMat$OR, na.rm = TRUE )
    tmpOR[ statMat$OR == -Inf ] <- min( statMat$OR, na.rm = TRUE )

    statMat$distance <- sapply(
        seq_along( statMat$Accession ),
        function(i){
            if( statMat$OR[i] > 1 ){
                x1 <- c( 0, 1 )
                x2 <- c( statMat$log10.adj.pVal[ i ], tmpOR[ i ] )
                d <- sqrt( sum( (x2-x1)**2 ) )
                
            } else if ( statMat$OR[ i ] <= 1 & statMat$OR[ i ] > 0 ) {
                x1 <- c( 0, 1 )
                x2 <- c( statMat$log10.adj.pVal[ i ], 1 / tmpOR[ i ] )
                d <- -1 * sqrt( sum( (x2-x1)**2 ) )
            } else { d <- 0 }
            
            return(d)
        }
    )

    if (!exists("MetaData")) {
        MetaData <- NULL
        data("MetaData", package = "TFEA.ChIP", envir = environment())
    }
    statMat <- merge(
        MetaData[, c( "Accession", "Cell", "Treatment" ) ],
        statMat, by = "Accession" )
    statMat <- statMat[ order( statMat$distance, decreasing = TRUE, na.last = TRUE ), ]
    return( statMat )
}

rankTFs <- function( resultsTable,
                     rankMethod = "gsea", makePlot = FALSE,
                     plotTitle = "TF ranking"){

    #' @title Rank the TFs in the output from 'getCMstats'
    #' @description Rank the TFs in the output from 'getCMstats' using
    #' Wilcoxon rank-sum test or a GSEA-like approach.
    #' @param resultsTable Output from the function 'getCMstats'
    #' @param rankMethod "wilcoxon" or "gsea".
    #' @param makePlot (Optional) For rankMethod="gsea". If TRUE, generates a plot for
    #' TFs with a p-value < 0.05.
    #' @param plotTitle (Optional) Title for the plot.
    #' @return data frame containing:
    #' \itemize{
    #'   \item For Wilcoxon rank-sum test: rank, TF name, test statistic
    #'   ('wilc_W), p-value, Freeman's theta, epsilon-squared anf effect size
    #'   \item For GSEA-like ranking: TF name, enrichment score, argument,
    #'   p-value, number of ChIPs}
    #' @export rankTFs
    #' @examples
    #' data('Genes.Upreg',package = 'TFEA.ChIP')
    #' CM_list_UP <- contingency_matrix(Genes.Upreg)
    #' stats_mat_UP <- getCMstats(CM_list_UP)
    #' rankTFs( stats_mat_UP )


    #### Input format
    rankMethod <- tolower( rankMethod )
    stopifnot(rankMethod %in% c("wilcoxon", "gsea"))

    # check the input type ( TFEA.ChIP association analysis  )
    cols_needed <- c("Accession","TF","distance")
    if( any( ! cols_needed %in% colnames( resultsTable ) ) ){
        stop("Input error: resultsTable doesn't contain all",
             "the columns required ('Accession','TF','distance')")
    }

    #### gsea ranking
    if( rankMethod == "gsea" ){

        chipSets <- lapply(
            unique( resultsTable$TF ),
            function(i, tab ){
                return( tab$Accession[ tab$TF == i ])
            }, tab = resultsTable
        )
        names( chipSets ) <- unique( resultsTable$TF )

        TFrank <- lapply(
            chipSets,
            function(i, statMat ){

                tf <-  statMat$TF[ statMat$Accession == i[1] ]
                res <- GSEA_EnrichmentScore(
                    statMat$Accession, i, 1, statMat$distance )

                shuffled <- rep( list( statMat$Accession ), 100)
                shuffled <- lapply( shuffled, sample )

                shuffled.ES <-  sapply(
                    shuffled,
                    function( j, i, chipDist ) {
                        tmp.ES <- GSEA_EnrichmentScore( j , i, 1, chipDist )$ES
                        return(tmp.ES)
                    }, i=i, chipDist = statMat$distance )

                shuffled.ES <- unlist( shuffled.ES )
                shuffled.ES <- shuffled.ES[ !is.na(shuffled.ES) ]
                pVal <- sum( abs(shuffled.ES) > abs(res$ES) ) / length( shuffled )

                return( data.frame(
                    "TF" = tf,
                    "ES" = res$ES,
                    "arg.ES" = res$arg.ES,
                    "pVal" = pVal,
                    "numberOfChIPs" = length(i),
                    stringsAsFactors = FALSE
                ))
            }, statMat = resultsTable
        )
        TFrank <- do.call( rbind, TFrank )

        if( makePlot == TRUE ){
            plot_df <- TFrank[TFrank$pVal<0.05 & !is.na(TFrank$pVal),]
            sub1 <- subset( plot_df, plot_df$ES > 0)
            sub2 <- subset( plot_df, plot_df$ES < 0)

            if( any( resultsTable$distance == 0 ) ){
                mid <- max(which( resultsTable$distance ==0 )) -
                    min( which( resultsTable$distance ==0 ) )/2
            } else {
                mid <- sum( resultsTable$distance > 0 ) + 0.5
            }

            p <- ggplot2::ggplot( ) +
              ggplot2::geom_point(
                  ggplot2::aes( x = arg.ES, y = ES, color = TF),
                  data = plot_df ) +
              ggplot2::ylim( -1.5, 1.5 ) +
              ggrepel::geom_text_repel( 
                  ggplot2::aes( x = arg.ES, y = ES, label = TF,
                                  color = TF),
                  data = sub1, nudge_y = 1,
                  angle = 90, direction = "x" ) +
              ggrepel::geom_text_repel(
                  ggplot2::aes( x = arg.ES, y = ES, label = TF,
                                  color = TF),
                  data = sub2, nudge_y = -1,
                  angle = 90, direction = "x" ) +
              ggplot2::theme_minimal() +
              ggplot2::geom_hline( yintercept = 0 ) +
              ggplot2::geom_point( ggplot2::aes( x=mid, y=0 ), 
                                   color="black" ) +
              ggplot2::guides( color = FALSE ) +
              ggplot2::labs( title = plotTitle,
                  y = "Enrichment Score", x= "ChIP-seq ranking")
            p
            return( list( TF_ranking=TFrank, TFranking_plot = p ))
        } else{ return(TFrank) }

        #### Wilcoxon rank-sum test
    } else if( rankMethod == "wilcoxon" ){

        chip_ranking <- data.frame(
            acc = resultsTable$Accession,
            TF = resultsTable$TF,
            rank = seq_along( resultsTable$Accession  ),
            stringsAsFactors = FALSE
        )

        TF_wilcox <- lapply(
            unique( chip_ranking$TF ),
            function(i, chip_ranking){

                tmp <- chip_ranking[, c(2:3)]
                tmp$TF[ tmp$TF != i ] <- "other"
                tmp$TF <- as.factor(tmp$TF)

                wilTest <- stats::wilcox.test( rank ~ TF, data=tmp )

                return( data.frame(
                    TF = i,
                    wilc_W = wilTest[["statistic"]],
                    wilc_pval = wilTest[["p.value"]],
                    freeman_theta = rcompanion::freemanTheta(x = tmp$rank, g = tmp$TF ),
                    epsilon_squared = rcompanion::epsilonSquared( x = tmp$rank, g = tmp$TF ),
                    wilc_r = rcompanion::wilcoxonR( x = tmp$rank, g = tmp$TF ),
                    stringsAsFactors = FALSE
                ))
            }, chip_ranking = chip_ranking
        )
        TF_wilcox <- do.call( rbind, TF_wilcox)
        TF_wilcox$rank <- seq_len( nrow( TF_wilcox ) )
        rownames( TF_wilcox ) <- TF_wilcox$rank
        TF_wilcox <- TF_wilcox[, c(7,1:6)]

        return( TF_wilcox )
    }
}

GSEA_EnrichmentScore <- function(gene.list, gene.set,
                                 weighted.score.type = 0, correl.vector = NULL) {

    # Computes the weighted GSEA score of gene.set in gene.list. Developed by
    # The Broad Institute

    #' @title Computes the weighted GSEA score of gene.set in gene.list.
    #' @description Computes the weighted GSEA score of gene.set in gene.list.
    #' @param gene.list The ordered gene list
    #' @param gene.set A gene set, e.g. gene IDs corresponding to a ChIP-Seq
    #' experiment's peaks.
    #' @param weighted.score.type Type of score: weight:
    #' 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
    #' @param correl.vector A vector with the coorelations (such as signal to
    #' noise scores) corresponding to the genes in the gene list
    #' @return list of:
    #' ES: Enrichment score (real number between -1 and +1)
    #' arg.ES: Location in gene.list where the peak running enrichment occurs
    #' (peak of the 'mountain')
    #' RES: Numerical vector containing the running enrichment score for all
    #' locations in the gene list
    #' tag.indicator: Binary vector indicating the location of the gene sets
    #' (1's) in the gene list
    #' @export GSEA_EnrichmentScore
    #' @examples
    #' GSEA_EnrichmentScore(gene.list=c('3091','2034','405','55818'),
    #' gene.set=c('2034','112399','405'))

    tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))
    no.tag.indicator <- 1 - tag.indicator
    N <- length(gene.list)
    Nh <- sum( gene.list %in% gene.set )
    Nm <- N - Nh
    if (weighted.score.type == 0) {
        correl.vector <- rep(1, N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector^alpha)
    sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
    if (sum.correl.tag > 0) {
        norm.tag <- 1/sum.correl.tag
        norm.no.tag <- 1/Nm
        RES <- cumsum(tag.indicator * correl.vector * norm.tag -
            no.tag.indicator * norm.no.tag)
        max.ES <- max( RES, na.rm = TRUE )
        min.ES <- min( RES, na.rm = TRUE )
        if (abs(max.ES) > abs(min.ES)) {
            # ES <- max.ES
            ES <- signif(max.ES, digits = 5)
            arg.ES <- which.max(RES)
        } else {
            # ES <- min.ES
            ES <- signif(min.ES, digits = 5)
            arg.ES <- which.min(RES)
        }
        return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))

    } else {
        RES <- rep(NaN, length(gene.list))
        indicator <- rep(0, length(gene.list))
        ES <- NA
        arg.ES <- NA
        return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
    }
}

GSEA_ESpermutations <- function(gene.list, gene.set, weighted.score.type = 0,
                                correl.vector = NULL, perms = 1000 ) {

    #' @title Calculate enrichment scores for a permutation test.
    #' @description Function to calculate enrichment scores over a randomly ordered
    #' gene list.
    #' @param gene.list Vector of gene Entrez IDs.
    #' @param gene.set A gene set, e.g. gene IDs corresponding to a ChIP-Seq
    #' experiment's peaks.
    #' @param weighted.score.type Type of score: weight:
    #' 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
    #' @param correl.vector A vector with the coorelations (such as signal to
    #' noise scores) corresponding to the genes in the gene list
    #' @param perms Number of permutations
    #' @return Vector of Enrichment Scores for a permutation test.
    # examples
    # GSEA_EnrichmentScore(gene.list=c('3091','2034','405','55818'),
    #' gene.set=c('2034','112399','405'), perms=10)

    tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))
    no.tag.indicator <- 1 - tag.indicator
    N <- length(gene.list)
    Nh <- sum( gene.list %in% gene.set )
    Nm <- N - Nh
    if (weighted.score.type == 0) {
        correl.vector <- rep(1, N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector^alpha)
    sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
    if (sum.correl.tag > 0) {
        norm.tag <- 1/sum.correl.tag
        norm.no.tag <- 1/Nm

        pES <- sapply(
            seq_len(perms),
            function( i, mInd, cv, mNorm, nmNorm){

                mInd <- sample( mInd )
                nmInd <- 1 - mInd
                sum.correl.tag <- sum(cv[mInd == 1])
                mNorm <- 1/sum.correl.tag

                RES <- cumsum(
                    mInd * cv * mNorm - nmInd * nmNorm)
                return( max( abs( RES ), na.rm = TRUE ) )
            },
            mInd=tag.indicator, cv=correl.vector,
            nmNorm=norm.no.tag
        )
        return( pES )

    } else {
        return( rep( NA, perms ) )
    }
}

GSEA_run <- function(gene.list, LFC, chip_index = get_chip_index(),
                     get.RES = FALSE, RES.filter = NULL, perms = 1000) {
  
    #' @title Function to run a GSEA analysis
    #' @description Function to run a GSEA to analyze the distribution of TFBS
    #' across a sorted list of genes.
    #' @param gene.list List of Entrez IDs ordered by their fold change.
    #' @param LFC Vector of log2( Fold Change ) values.
    #' @param chip_index Output of the function “get_chip_index”, a data frame
    #' containing accession IDs of ChIPs on the database and the TF each one
    #' tests. If not provided, the whole internal database will be used
    #' @param get.RES (Optional) boolean. If TRUE, the function stores Running
    #' Enrichment Scores of all/some TF.
    #' @param RES.filter (Optional) chr vector. When get.RES==TRUE, allows to
    #' choose which TF's Running Enrichment Score to store.
    #' @param perms (Optional) integer. Number of permutations for the enrichment test.
    #' @return a list of:
    #' Enrichment.table: data frame containing accession ID, Cell type, ChIP-Seq
    #' treatment, transcription factor tested, enrichment score, adjusted p-value,
    #' and argument of every ChIP-Seq experiment.
    #' RES (optional): list of running sums of every ChIP-Seq
    #' indicators (optional): list of 0/1 vectors that stores the matches (1)
    #' and mismatches (0) between the gene list and the gene set.
    #' @export GSEA_run
    #' @examples
    #' data( 'hypoxia', package = 'TFEA.ChIP' )
    #' hypoxia <- preprocessInputData( hypoxia )
    #' chip_index <- get_chip_index( TFfilter = c('HIF1A','EPAS1','ARNT' ) )
    #' GSEA.result <- GSEA_run( hypoxia$Genes, hypoxia$log2FoldChange, chip_index, get.RES = TRUE)
  
    if (!exists("ChIPDB")) {
        ChIPDB <- NULL
        data("ChIPDB", package = "TFEA.ChIP", envir = environment())
    }
    if( is.matrix( ChIPDB ) ){
        ChIPDB <- matrixDB_to_listDB( ChIPDB )
    }
    if (!exists("MetaData")) {
        MetaData <- NULL
        data("MetaData", package = "TFEA.ChIP", envir = environment())
    }
    
    enrichmentScore <- vector()
    pval <- vector()
    enrichmentArg <- vector()
    
    if (get.RES == TRUE) {
        res <- list()
        ind <- list()
    }
    # create progress bar
    pbar <- txtProgressBar(
        min = 0, max = length(chip_index$Accession), style = 3)
    
    for (i in seq_along(chip_index$Accession)) {
        
        acc <- chip_index$Accession[i] 
        targets <- ChIPDB[["Gene Keys"]][ ChIPDB[["ChIP Targets"]][[acc]] ]
        
        if ( sum( targets %in% gene.list ) > 10 ) {
            
            result <- GSEA_EnrichmentScore( gene.list, targets, 1, LFC )
            
            shuffled.ES <- GSEA_ESpermutations( 
                gene.list , targets, 1, LFC, perms )
            shuffled.ES <- shuffled.ES[ !is.na( shuffled.ES ) ]
            
            enrichmentScore <- c( enrichmentScore, result$ES )
            pval <- c( pval,
                       sum( shuffled.ES >= abs(result$ES))/ length(shuffled.ES))
            enrichmentArg <- c( enrichmentArg, result$arg.ES )
            
            if( get.RES ==TRUE ){
                
                if ( missing( RES.filter ) ) {
                    # Store running sums of all TFs.
                    res <- c( res, list( result$RES ) )
                    names( res )[ length( res )] <- acc
                    ind <- c( ind, list( result$indicator ) )
                    names( ind )[ length( ind ) ] <- acc
                    
                } else {
                    # Store running sums of TF in RES.filter
                    if ( acc %in% RES.filter) {
                        res <- c( res, list( result$RES ) )
                        names( res )[ length( res )] <- acc
                        ind <- c( ind, list( result$indicator ) )
                        names( ind )[ length( ind )] <- acc
                    }
                }
            }
            
        } else {
            chip_index <- chip_index[-i, ]
            i <- i-1
        }
        # update progress bar
        setTxtProgressBar(pbar, i)
    }
    close(pbar)
    pval.adj <- stats::p.adjust(pval, "fdr")  # Adjust pvalues
    
    enrichmentTable <- data.frame(
        Accession = chip_index$Accession, TF = chip_index$TF,
        ES = enrichmentScore, p.val = pval, pval.adj = pval.adj,
        Arg.ES = enrichmentArg, stringsAsFactors = F )
    
    enrichmentTable <- merge(
        MetaData[,c("Accession","Cell","Treatment")],
        enrichmentTable, by="Accession")
    
    enrichmentTable <- enrichmentTable[ !is.na( enrichmentTable$pval.adj ), ]
    
    if ( get.RES == TRUE) {
        GSEA_results <- list( enrichmentTable, res, ind )
        names( GSEA_results ) <- c("Enrichment.table", "RES", "indicators")
        return( GSEA_results )
    } else {
        return( enrichmentTable )
    }
}

#### Plotting functions ####

plot_CM <- function( CM.statMatrix, plot_title = NULL,
    specialTF = NULL, TF_colors = NULL ) {
    
    #' @title Makes an interactive html plot from an enrichment table.
    #' @description Function to generate an interactive html plot from a
    #' transcription factor enrichment table, output of the function
    #' 'getCMstats'.
    #' @param CM.statMatrix Output of the function 'getCMstats'.
    #' A data frame storing: Accession ID of every ChIP-Seq tested,
    #' Transcription Factor,Odds Ratio, p-value and adjusted p-value.
    #' @param plot_title The title for the plot.
    #' @param specialTF (Optional) Named vector of TF symbols -as written in
    #' the enrichment table- to be highlighted in the plot. The name of each
    #' element of the vector specifies its color group, i.e.: naming elements
    #' HIF1A and HIF1B as 'HIF' to represent them with the same color.
    #' @param TF_colors (Optional) Nolors to highlight TFs chosen in specialTF.
    #' @return plotly scatter plot.
    #' @export plot_CM
    #' @examples
    #' data('Genes.Upreg',package = 'TFEA.ChIP')
    #' CM_list_UP <- contingency_matrix( Genes.Upreg )
    #' stats_mat_UP <- getCMstats( CM_list_UP )
    #' plot_CM( stats_mat_UP )
    
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("plotly package needed for this function to work. ",
             "Please install it.", call. = FALSE)
    }
    
    # Checking input variables
    if (is.null(plot_title)) {
        plot_title <- "Transcription Factor Enrichment"
    }
    if( is.null( specialTF ) ) {
        CM.statMatrix$highlight <- rep("Other", dim(CM.statMatrix)[1])
        markerColors <- c("azure4")
        names(markerColors) <- c("Other")
    }
    if( !is.null( specialTF ) & is.null( TF_colors ) ) {
        TF_colors <- c(
            "red", "blue", "green", "hotpink", "cyan", "greenyellow", "gold",
            "darkorchid", "chocolate1", "black", "lightpink", "seagreen")
        TF_colors <- TF_colors[ 1 : length( unique( names( specialTF ) ) ) ]
        color_list <- highlight_TF( CM.statMatrix, 4, specialTF, TF_colors ) 
        CM.statMatrix$highlight <- color_list[[ 1 ]]
        markerColors <- color_list[[ 2 ]]
    }
    if (!is.null(specialTF) & !is.null(TF_colors)) {
        color_list <- highlight_TF( CM.statMatrix, 4, specialTF, TF_colors )
        CM.statMatrix$highlight <- color_list[[ 1 ]]
        markerColors <- color_list[[ 2 ]]
    }
    if (!exists("MetaData")) {
        MetaData <- NULL
        data("MetaData", package = "TFEA.ChIP", envir = environment())
    }
    
    # Adding metadata for the plot
    MetaData <- MetaData[
        match( CM.statMatrix$Accession, MetaData$Accession), ]
    CM.statMatrix$Treatment <- MetaData$Treatment
    CM.statMatrix$Cell <- MetaData$Cell
    rm(MetaData)
    
    # Cheking if any plot variables have a Inf values.
    if ( sum( CM.statMatrix$OR == Inf ) > 0 ) {
        
        warn_number <- sum( CM.statMatrix$OR == Inf )
        
        # Substitute Inf OR values for the maximum finite value
        finMax <- max( CM.statMatrix$OR[ CM.statMatrix$OR != Inf ] )
        CM.statMatrix[CM.statMatrix$OR == Inf, ]$OR <- finMax
        
        warning(warn_number, " elements have an Odds Ratio of Inf.",
                " Maximum value for OR introduced instead.")
    }
    if ( sum( CM.statMatrix$OR == -Inf ) > 0) {
        
        warn_number <- sum( CM.statMatrix$OR == -Inf )
        
        # Substitute -Inf OR values for the minimum finite value
        finMin <- min( CM.statMatrix$OR[ CM.statMatrix$OR != -Inf ] )
        CM.statMatrix$OR[CM.statMatrix$OR == -Inf, ]$OR <- finMin
        
        warning(warn_number, " elements have an Odds Ratio of -Inf. Minimum",
                " value for OR introduced instead.")
    }
    if ( sum( CM.statMatrix$adj.p.value == 0 ) > 0) {
        warn_number <- sum( CM.statMatrix$adj.p.value == 0 )
        
        # Substitute Inf -log(pval) values for the maximum finite value
        finMax <- max( CM.statMatrix$log.adj.pVal[ CM.statMatrix$p.value != 0])
        CM.statMatrix$log.adj.pVal[ CM.statMatrix$p.value == 0] <- finMax
        
        warning(warn_number, " elements have a -log(p-Value) of Inf. ",
                "Maximum value for -log(p-Val) introduced instead.")
    }
    
    CM.statMatrix$pointText <- paste0(
        CM.statMatrix$Accession, ": ", CM.statMatrix$TF,
        "<br>Treatment: ", CM.statMatrix$Treatment,
        "<br>Cell: ", CM.statMatrix$Cell)
    
    if (length(markerColors) > 1) {
        colorPoints <- CM.statMatrix[ CM.statMatrix$highlight != "Other", ]
        bgPoints <- CM.statMatrix[ CM.statMatrix$highlight == "Other", ]
        
        p <- plotly::plot_ly( bgPoints, 
            x = ~log10.adj.pVal, y = ~log2.OR, type = "scatter",
            mode = "markers", text = ~pointText, color = ~highlight,
            colors = markerColors )
        
        p <- plotly::add_markers( p,
            x = colorPoints$log10.adj.pVal, y = colorPoints$log2.OR,
            type = "scatter", mode = "markers", text = colorPoints$pointText,
            color = colorPoints$highlight, colors = markerColors ) %>%
            plotly::layout( title = plot_title )
        
    } else {
        p <- plotly::plot_ly( CM.statMatrix, 
            x = ~log10.adj.pVal, y = ~log2.OR, type = "scatter", 
            mode = "markers", text = ~pointText, color = ~highlight,
            colors = markerColors ) %>%
            plotly::layout( title = plot_title )
    }
    p
    return(p)
}

plot_ES <- function( GSEA_result, LFC, plot_title = NULL, specialTF = NULL,
    TF_colors = NULL, Accession = NULL, TF = NULL ) {
    
    #' @title Plots Enrichment Score from the output of GSEA.run.
    #' @description Function to plot the Enrichment Score of every member of
    #' the ChIPseq binding database.
    #' @param GSEA_result Returned by GSEA_run
    #' @param LFC Vector with log2(Fold Change) of every gene that has an
    #' Entrez ID. Arranged from higher to lower.
    #' @param plot_title (Optional) Title for the plot
    #' @param specialTF (Optional) Named vector of TF symbols -as written in
    #' the enrichment table- to be highlighted in the plot. The name of each
    #' element specifies its color group, i.e.: naming elements HIF1A and HIF1B
    #' as 'HIF' to represent them with the same color.
    #' @param TF_colors (Optional) Colors to highlight TFs chosen in specialTF.
    #' @param Accession (Optional) restricts plot to the indicated list dataset
    #' IDs.
    #' @param TF (Optional) restricts plot to the indicated list transcription
    #' factor names.
    #' @return Plotly object with a scatter plot -Enrichment scores- and a
    #' heatmap -log2(fold change) bar-.
    #' @export plot_ES
    #' @examples
    #' data('GSEA.result','log2.FC',package = 'TFEA.ChIP')
    #' TF.hightlight <- c('E2F1' = 'E2F1')
    #' col <- c('red')
    #' plot_ES( GSEA.result, log2.FC, "Example", TF.hightlight, col )
    
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("plotly package needed for this function to work. ",
             "Please install it.", call. = FALSE)
    }
    
    if ( is.data.frame( GSEA_result ) ) {
        enrichTab <- GSEA_result
    } else {
        enrichTab <- GSEA_result[["Enrichment.table"]]
    }
    
    # If plotting selected points only
    if ( ! is.null( Accession ) | ! is.null( TF ) ) {
        if ( is.null( Accession ) ) {
            Accession <- enrichTab$Accession[ enrichTab$TF %in% TF ]
        }
        if (is.null(TF)) {
            TF <- enrichTab$TF[ enrichTab$Accession %in% Accession ]
        }
        sel <- enrichTab$Accession %in% Accession & enrichTab$TF %in% TF
        enrichTab <- enrichTab[ which( sel ), ]
    }
    
    # Default options
    if ( is.null(plot_title) ) {
        plot_title <- "Transcription Factor Enrichment"
    }
    
    if ( is.null( specialTF ) ) {
        enrichTab$highlight <- "Other"
        markerColors <- c( "Other"="azure4")
    } else {
        if ( is.null( TF_colors ) ) {
            TF_colors <- c(
                "red", "blue", "green", "hotpink", "cyan", "greenyellow",
                "gold", "darkorchid", "chocolate1", "black", "lightpink",
                "seagreen" )
            TF_colors <- TF_colors[ 1 : length( unique( names( specialTF ) ) ) ]
            color_list <- highlight_TF( enrichTab, 4, specialTF, TF_colors )
            enrichTab$highlight <- color_list[[ 1 ]]
            markerColors <- color_list[[ 2 ]]
        } else {
            color_list <- highlight_TF( enrichTab, 4, specialTF, TF_colors )
            enrichTab$highlight <- color_list[[ 1 ]]
            markerColors <- color_list[[ 2 ]]
        }
    }

    enrichTab$symbol <- "pVal>0.05"
    enrichTab$symbol[ enrichTab$pval.adj <= 0.05 ] <- "pVal<=0.05"
    
    # Adding metadata
    if ( ! exists( "MetaData" ) ) {
        MetaData <- NULL
        data("MetaData", package = "TFEA.ChIP", envir = environment())
    }
    
    MetaData <- MetaData[ match( enrichTab$Accession, MetaData$Accession), ]
    enrichTab$Treatment <- MetaData$Treatment
    enrichTab$Cell <- MetaData$Cell
    rm(MetaData)
    
   
    enrichTab$pointText = paste0(
        enrichTab$Accession, ": ", enrichTab$TF, "<br>Adjusted p-value: ",
        round( enrichTab$pval.adj, 3 ), "<br>Treatment: ", enrichTab$Treatment,
        "<br>Cell: ", enrichTab$Cell )
    
    
    multicolor <- length(markerColors) > 1 &
        length( unique( enrichTab$TF ) ) > length( specialTF )
    
    if ( multicolor ) {
        colorPoints <- enrichTab[ enrichTab$highlight != "Other", ]
        bgPoints <- enrichTab[ enrichTab$highlight == "Other", ]
        
        p <- plotly::plot_ly( bgPoints,
            x = bgPoints$Arg.ES, y = bgPoints$ES, type = "scatter",
            mode = "markers", text = bgPoints$pointText,
            color = bgPoints$highlight, colors = markerColors,
            symbol = bgPoints$symbol, symbols = c( "x", "circle" ) )
        
        p <- plotly::add_markers( p, 
            x = colorPoints$Arg.ES, y = colorPoints$ES, type = "scatter",
            mode = "markers", text = colorPoints$pointText,
            color = colorPoints$highlight, colors = markerColors,
            symbol = colorPoints$symbol, symbols = c("x", "circle") ) %>%
            plotly::layout( title = plot_title,
                xaxis = list(title = "Argument"), yaxis = list(title = "ES") )
        
    } else {
        p <- plotly::plot_ly(enrichTab,
            x = enrichTab$Arg.ES, y = enrichTab$ES, type = "scatter", 
            mode = "markers", text = enrichTab$pointText, 
            color = enrichTab$highlight, colors = markerColors, 
            symbol = enrichTab$symbol, symbols = c( "x", "circle" ) ) %>%
            plotly::layout( title = plot_title,
                xaxis = list( title = "Argument"), yaxis = list(title = "ES"))
    }
    # Adding log2(Fold Change) bar to the plot
    LFCbar <- get_LFC_bar( LFC )
    graf <- plotly::subplot( p, 
        LFCbar, shareX = TRUE, nrows = 2, heights = c(0.95, 0.05),
        titleY = TRUE)
    graf
    return(graf)
}

plot_RES <- function( GSEA_result, LFC, plot_title = NULL, line.colors = NULL,
    line.styles = NULL, Accession = NULL, TF = NULL ) {

    #' @title Plots all the RES stored in a GSEA_run output.
    #' @description Function to plot all the RES stored in a GSEA_run output.
    #' @param GSEA_result Returned by GSEA_run
    #' @param LFC Vector with log2(Fold Change) of every gene that has an
    #' Entrez ID. Arranged from higher to lower.
    #' @param plot_title (Optional) Title for the plot.
    #' @param line.colors (Optional) Vector of colors for each line.
    #' @param line.styles (Optional) Vector of line styles for each line
    #' ('solid'/'dash'/'longdash').
    #' @param Accession (Optional) restricts plot to the indicated list dataset
    #' IDs.
    #' @param TF (Optional) restricts plot to the indicated list transcription
    #' factor names.
    #' @return Plotly object with a line plot -running sums- and a
    #' heatmap -log2(fold change) bar-.
    #' @export plot_RES
    #' @examples
    #' data('GSEA.result','log2.FC',package = 'TFEA.ChIP')
    #' plot_RES(GSEA.result, log2.FC, TF = c('E2F4',"E2F1"),
    #'     Accession=c('ENCSR000DYY.E2F4.GM12878',
    #'     'ENCSR000EVJ.E2F1.HeLa-S3'))

    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("plotly package needed for this function to work.",
            " Please install it.", call. = FALSE)
    }
    # Checking input variables
    if (!exists("MetaData")) {
        MetaData <- NULL
        data("MetaData", package = "TFEA.ChIP", envir = environment())
    }
    
    enrichTab <- GSEA_result[["Enrichment.table"]]
    rSums <- GSEA_result[["RES"]]
    rm( GSEA_result )
    
    # Accession or TF filters
    if ( !is.null( Accession ) | !is.null( TF ) ) {

        if( is.null( Accession ) ) { 
            Accession <- enrichTab$Accession[ enrichTab$TF %in% TF ] }
        if( is.null( TF ) ) { 
            TF <- enrichTab$TF[ enrichTab$Accession %in% Accession ] }

        enrichTab <- enrichTab[ match( Accession, enrichTab$Accession ), ]
        rSums <- rSums[ match( Accession, names( rSums ) ) ]

    } else { Accession <- enrichTab$Accession }
    
    # line parameters & plot title
    if( is.null( line.colors ) ) {
        line.colors <- c("red", "blue", "green", "hotpink", "cyan",
            "greenyellow", "gold", "darkorchid", "chocolate1",
            "black", "lightpink", "seagreen")
        line.colors <- line.colors[ 1 : length( Accession ) ]
    }
    if( is.null( line.styles ) ) {
        line.styles <- rep("solid", length( Accession ))
    }
    if( is.null( plot_title ) ) {
        plot_title <- "Transcription Factor Enrichment"
    }

    # Getting metadata for the plot
    MetaData <- MetaData[ match( Accession, MetaData$Accession ), ]

    plotTab <- data.frame(
    	Accession = MetaData$Accession,
    	Cell = MetaData$Cell,
    	Treatment = MetaData$Treatment,
    	TF = MetaData$TF,
    	stringsAsFactors = FALSE)

    plotTab$RES <- rSums
    rm( MetaData, rSums )
    
    plotTab$lineName <- paste0( plotTab$Accession, " - ", plotTab$TF )
    plotTab$lineText <- paste0(
        plotTab$Accession, " - ", plotTab$TF,
        "<br>Cell: ", plotTab$Cell, "<br>Treatment: ", plotTab$Treatment )

    if ( dim( plotTab )[ 1 ] > 1 ) { # plot multiple lines
        for (i in seq_along(Accession)) {
            if (i == 1) { # First line
                grafica <- plotly::plot_ly( plotTab,
                    x = c( 1 : length( plotTab$RES[[ 1 ]] ) ),
                    y = plotTab$RES[[ Accession[ 1 ] ]], 
                    type = "scatter", mode = "lines", 
                    line = list(color = line.colors[1], dash = line.styles[1]),
                    name = plotTab$lineName[ 1 ], text = plotTab$lineText[ 1 ])
                
            } else if ( i > 1 & i < length( Accession ) ) {
                grafica <- plotly::add_trace( grafica,
                    y = plotTab$RES[[ Accession[ i ] ]],
                    type = "scatter", mode = "lines",
                    line = list(color = line.colors[i], dash = line.styles[i]),
                    name = plotTab$lineName[ i ], text = plotTab$lineText[ i ])
                
            } else if (i == length(Accession)) { # Last line
                grafica <- plotly::add_trace( grafica,
                    y = plotTab$RES[[Accession[i]]],
                    type = "scatter", mode = "lines",
                    line = list(color = line.colors[i], dash = line.styles[i]),
                    name = plotTab$lineName[ i ],
                    text = plotTab$lineText[ i ] ) %>%
                plotly::layout(
                    title = plot_title, xaxis = list(title = "Argument"),
                    yaxis = list(title = "ES"))
            }
        }
    } else { # plot one line
        grafica <- plotly::plot_ly( plotTab, 
            x = c( 1 : length( plotTab$RES[[ 1 ]] ) ),
            y = plotTab$RES[[ Accession[ 1 ] ]], 
            type = "scatter", mode = "lines",
            line = list(color = line.colors[1], dash = line.styles[1]),
            name = plotTab$lineText ) %>%
        plotly::layout(title = plot_title,
            xaxis = list(title = "Argument"), yaxis = list(title = "ES"))
    }

    # Adding log2(Fold Change) bar to the plot
    LFC.bar <- get_LFC_bar( LFC )
    graf <- plotly::subplot( grafica, LFC.bar, shareX = TRUE,
        nrows = 2, heights = c(0.95, 0.05), titleY = TRUE )
    graf
    return( graf )
}

highlight_TF <- function(table, column, specialTF, markerColors) {

    #' @title Highlight certain transcription factors in a plotly graph.
    #' @description Function to highlight certain transcription factors using
    #' different colors in a plotly graph.
    #' @param table Enrichment matrix/data.frame.
    #' @param column Column # that stores the TF name in the matrix/df.
    #' @param specialTF Named vector containing TF names as they appear in the
    #' enrichment matrix/df and nicknames for their color group.
    #' Example:
    #'           specialTF<-c('HIF1A','EPAS1','ARNT','SIN3A')
    #'           names(specialTF)<-c('HIF','HIF','HIF','SIN3A')
    #' @param markerColors Vector specifying the shade for every color group.
    #' @return List of two objects:
    #' A vector to attach to the enrichment matrix/df pointing out the color
    #' group of every row.
    #' A named vector connecting each color group to the chosen color.
    # examples
    # highlight_TF( CM.statMatrix_UP, 4, specialTF, colors )

    highlight <- sapply( seq_along(table[,1]),
        function (table, column, TF, i) {
            if (table[i, column] %in% TF){
                color.group <- names( TF[ TF == table[i, column] ] )
            }else{
                color.group <- "Other"
            }
            return(color.group) },
        table = table, column = column, TF = specialTF
        )

    markerColors <- c( "azure4", markerColors )
    names( markerColors ) <- c( "Other", unique( names( specialTF ) ) )
    return( list( highlight, markerColors ) )
}

get_LFC_bar <- function( LFC ) {

    #' @title Plots a color bar from log2(Fold Change) values.
    #' @description Function to plot a color bar from log2(Fold Change)
    #' values from an expression experiment.
    #' @param LFC Vector of log2(fold change) values arranged from higher
    #' to lower. Use ony the values of genes that have an Entrez ID.
    #' @return Plotly heatmap plot -log2(fold change) bar-.
    # examples
    # get_LFC_bar( arranged.log2FC.array )

    if (!requireNamespace("scales", quietly = TRUE)) {
        stop("scales package needed for this function to work. ",
            "Please install it.", call. = FALSE)
    }
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("plotly package needed for this function to work. ",
            "Please install it.", call. = FALSE)
    }

    # Rescaling log Fold Change values to get a -1 to 1 scale
    vals <- scales::rescale( LFC )
    o <- order( vals, decreasing = FALSE )

    posLen <- sum( LFC > 0 )

    # Building a red-blue scale adapted to the log Fold Change
    cols1 <- ( scales::col_numeric(
        grDevices::colorRamp( c( "mistyrose", "red3" ) ),
        domain = NULL) )( vals[ 1 : posLen ] )
    
    cols2 <- ( scales::col_numeric(
        grDevices::colorRamp(c("navy", "lightcyan")),
        domain = NULL)
        )( vals[ posLen+1 : length(LFC) ] )
    
    cols <- c(cols1, cols2)

    colorValues <- data.frame( vals[ o ], cols[ o ] )

    LFC.bar <- plotly::plot_ly(
        x = c( 1 : length( LFC ) ), y = rep( 1, length( LFC ) ), z = LFC,
        type = "heatmap", colorscale = colorValues, showscale = FALSE) %>%
    plotly::layout(yaxis = list(visible = FALSE))

    return( LFC.bar )
}
