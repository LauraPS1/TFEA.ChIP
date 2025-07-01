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
#' @import org.Mm.eg.db
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
      # adjust p-values using Benjamini & Hochberg correction
      fileTable$score <- p.adjust( fileTable$score, "fdr" ) 
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
                                   ranges = IRanges::IRanges(fileTable$start, 
                                                          end = fileTable$end),
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
        # -10 * log10(p-value)
        fileTable$score <- 10 ^ (-1 * (fileTable$score / 10 ) ) 
      } else if (format == "macs2"){
        fileTable$score <- 10 ^ (-1 * (fileTable$score ) ) # -log10(p-value)
      }
      # adjust p-values using Benjamini & Hochberg correction
      fileTable$score <- p.adjust( fileTable$score, "fdr" ) 
      
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
        # -10 * log10(p-value)
        fileTable$score <- 10 ^ (-1 * (fileTable$score / 10 ) ) 
      } else if (format == "macs2"){
        fileTable$score <- 10 ^ (-1 * (fileTable$score ) ) # -log10(p-value)
      }
      
      # adjust p-values using Benjamini & Hochberg correction
      fileTable$score <- p.adjust( fileTable$score, "fdr" ) 
      fileTable <- fileTable[fileTable$score < valLimit, ]
      Stat <- "corrected p-Value"
    }
    
    if( dim( fileTable)[1] > 0){
      
      fileMetaData <- c(fileMetaData, Stat)
      MDframe <- as.data.frame(lapply(fileMetaData, rep, dim(fileTable)[1]))
      colnames(MDframe) <- c("Name", "Accession", "Cell", "Cell Type",
                             "Treatment", "Antibody", "TF", "Score Type")
      
      gr <- GenomicRanges::GRanges(seqnames = fileTable$chr,
                                   ranges = IRanges::IRanges(fileTable$start, 
                                   end = fileTable$end),
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
  #' # For this example, we will use the variables already included in the
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
  #' is sorted according to decreasing log2(Fold Change). Translating
  #' gene IDs from mouse to their equivalent human genes is available
  #' using the variable "mode".
  #' @param inputData DESeqResults object or data frame. In all cases
  #' must include gene IDs. Data frame inputs should include 'pvalue' and
  #' 'log2FoldChange' as well.
  #' @param  mode Specify the organism used: 'h2h' for homo sapiens gene IDs,
  #' 'm2m' for mouse gene IDs, or 'm2h' to get the corresponding human gene
  #' IDs from a mouse input.
  #' @return A table containing Entrez Gene IDs, Gene Symbols, LogFoldChange 
  #' and p-val values (both raw p-value and fdr adjusted p-value), sorted by
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
        gene_symbols <- genes$GENE.ID
        entrez_ids <- genes$ENTREZ.ID
      } else {
        genes <- genes[!is.na(genes$human.gene.ID), ]
        inputData <- inputData[rownames(inputData) %in% genes$mouse.gene.ID, ]
        gene_symbols <- genes$mouse.gene.ID
        entrez_ids <- genes$human.gene.ID
      }
      
    } else {
      gene_symbols <- rownames(inputData)
      # No translation needed if already in Entrez format
      entrez_ids <- gene_symbols  
    }
    # get the rest of the variables
    log2FoldChange <- inputData[["log2FoldChange"]]
    pvalue <- inputData[["pvalue"]]
    pval.adj <- inputData[["padj"]]
    
    Table <- data.frame(Symbol = gene_symbols, Genes = entrez_ids, 
                        log2FoldChange = log2FoldChange, 
                        pvalue = pvalue, pval.adj = pval.adj)
    Table$Genes <- as.character(Table$Genes)
    Table <- Table[!is.na(Table$log2FoldChange), ]
    Table <- Table[order(Table$log2FoldChange, decreasing = TRUE), ]
    rownames(Table) <- NULL
    return(Table)
    
  } else if ( methods::is(inputData, "TopTags") ){
    # Extracting data from a TopTags object
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("edgeR package needed for this function to work. ",
           "Please install it.", call. = FALSE)
    }
    inputData <- as.data.frame( inputData )
    
    # check the gene ids and translate if needed
    if ( ! all( grepl("^\\d*$", rownames(inputData) ) ) | mode == "m2h" ) {
      genes <- suppressMessages( GeneID2entrez(
        gene.IDs = rownames( inputData ),
        mode,
        return.Matrix = TRUE))
      
      if ( mode %in% c("h2h","m2m") ){
        genes <- genes[!is.na(genes$ENTREZ.ID), ]
        inputData <- inputData[rownames(inputData) %in% genes$GENE.ID, ]
        gene_symbols <- genes$GENE.ID
        entrez_ids <- genes$ENTREZ.ID
      } else {
        genes <- genes[!is.na(genes$human.gene.ID), ]
        inputData <- inputData[rownames(inputData) %in% genes$mouse.gene.ID, ]
        gene_symbols <- genes$mouse.gene.ID
        entrez_ids <- genes$human.gene.ID
      }
      
    } else {
      gene_symbols <- rownames(inputData)
      # No translation needed if already in Entrez format
      entrez_ids <- gene_symbols  
    }
    # get the rest of the variables
    log2FoldChange <- inputData$logFC
    pvalue <- inputData$PValue
    pval.adj <- inputData$FDR
    
    Table <- data.frame(Symbol = gene_symbols, Genes = entrez_ids, 
                        log2FoldChange = log2FoldChange, pvalue = pvalue, 
                                                        pval.adj = pval.adj)
    Table$Genes <- as.character(Table$Genes)
    Table <- Table[!is.na(Table$log2FoldChange), ]
    Table <- Table[order(Table$log2FoldChange, decreasing = TRUE), ]
    rownames(Table) <- NULL
    return(Table)
    
  } else if ( methods::is(inputData, "data.frame") ) {
    # Extracting data from a data frame
    # Checkig if all the necessary columns are present
    if (!(all(c("Genes", "pvalue", "log2FoldChange") %in% colnames(inputData))) && 
        !(all(c("Genes", "pval.adj", "log2FoldChange") %in% colnames(inputData)))) {
      stop("Input data must include: 'Genes', 'log2FoldChange', and either 'pvalue' or 'pval.adj'.")
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
        inputData$Symbol <- genes$GENE.ID
        inputData$Genes <- genes$ENTREZ.ID
      }else{
        genes <- genes[!is.na(genes$human.gene.ID), ]
        inputData <- inputData[ inputData$Genes %in% genes$mouse.gene.ID, ]
        inputData$Symbol <- genes$mouse.gene.ID
        inputData$Genes <- genes$human.gene.ID
      }
    }
    
    # remove NA entrez rows
    inputData <- inputData[!is.na(inputData$Genes), ]
    
    # make sure log2(FoldChange) are numbers
    inputData$log2FoldChange <- as.numeric(inputData$log2FoldChange)
    
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
  #' @param GeneExpression_df A data frame with the following fields:
  #' 'Gene', 'pvalue' or 'pval.adj', 'log2FoldChange'.
  #' @param max_pval maximum p-value allowed, 0.05 by default.
  #' @param min_pval minimum p-value allowed, 0 by default.
  #' @param max_LFC maximum log2(FoldChange) allowed.
  #' @param min_LFC minimum log2(FoldChange) allowed.
  #' @return A character vector of gene IDs.
  #' @export Select_genes
  #' @examples
  #' data('hypoxia', package='TFEA.ChIP')
  #' Select_genes(hypoxia)
  
  # Checking input variables
  if (max_pval < min_pval) {
    stop("'max_pval' must be greater than or equal to 'min_pval'.",
                                                         call. = FALSE)
  }
  if (max_LFC < min_LFC) {
    stop("'max_LFC' must be greater than or equal to 'min_LFC'.", 
                                                         call. = FALSE)
  }
  
  # Ensure required columns are present
  if (!all(c("Genes", "log2FoldChange") %in% colnames(GeneExpression_df))) {
    stop("The input data must contain 'Genes' and 'log2FoldChange'.", 
                                                           call. = FALSE)
  }
  
  # Ensure at least one p-value column is present
  if (!any(c("pval.adj", "pvalue") %in% colnames(GeneExpression_df))) {
    stop("The input data must contain either 'pval.adj' or 'pvalue'.", 
                                                             call. = FALSE)
  }
  
  # Selecting genes based on p-value and log2(FoldChange)
  if ("pval.adj" %in% colnames(GeneExpression_df)) {
    g_pv <- with(GeneExpression_df, pval.adj >= min_pval & pval.adj <= max_pval)
  } else {
    g_pv <- with(GeneExpression_df, pvalue >= min_pval & pvalue <= max_pval)
  }
  
  g_lfc <- with(GeneExpression_df,
                      log2FoldChange >= min_LFC & log2FoldChange <= max_LFC)
  
  # Extracting gene names
  g_names <- GeneExpression_df$Genes[g_pv & g_lfc]
  
  return(as.character(g_names))
}


GeneID2entrez <- function(gene.IDs, return.Matrix = FALSE, mode = "h2h") {
  
  #' @title Translates gene IDs from Gene Symbol or Ensembl ID to Entrez ID.
  #' @description Translates mouse or human gene IDs from Gene Symbol or
  #' Ensembl Gene ID to Entrez Gene ID using AnnotationDbi.
  #' @param gene.IDs Array of Gene Symbols or Ensembl Gene IDs.
  #' @param return.Matrix Logical. When TRUE, the function returns a matrix [n,2],
  #' one column with the gene symbols or Ensembl IDs, another with their
  #' respective Entrez IDs.
  #' @param mode Specify the organism used: 'h2h' for human gene IDs,
  #' 'm2m' for mouse gene IDs, or 'm2h' for converting mouse to human gene IDs.
  #' @return Vector or matrix containing the Entrez IDs (or NA) corresponding
  #' to every element of the input.
  #' @export GeneID2entrez
  #' @examples
  #' GeneID2entrez(c('TNMD','DPM1','SCYL3','FGR','CFH','FUCA2','GCLC'))
  
  # Ensure valid mode
  stopifnot(mode %in% c("h2h", "m2m", "m2h"))
  
  gene.IDs <- gene.IDs[!is.na(gene.IDs)]
  gene.IDs <- trimws(gene.IDs)  # Remove any possible white spaces
  
  # Helper function to handle warnings for many-to-one mappings
  handle_warning <- function(matched_2) {
    if (sum(duplicated(matched_2[!is.na(matched_2)])) > 0) {
      warning("Some genes returned 1:many mapping to ENTREZ ID.",
              call. = FALSE)
    }
  }
  
  if (mode == 'h2h') {
    gene.IDs <- toupper(gene.IDs)  # Ensure uppercase for human symbols
    
    # Check if gene.IDs are Ensembl or Gene Symbols
    if (all(grepl("^ENSG", gene.IDs, perl = TRUE))) {
      ID.type <- "ENSEMBL"
      suppressMessages(GeneNames <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db, gene.IDs, c("ENTREZID", "ENSEMBL"), 
        keytype = "ENSEMBL"))
    } else {
      ID.type <- "SYMBOL"
      suppressMessages(GeneNames <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db, gene.IDs, c("SYMBOL", "ENTREZID"), 
        keytype = "ALIAS"))
    }
    
    # Match results
    matched <- match(as.character(gene.IDs), GeneNames[, ID.type])
    matched_2 <- match(GeneNames[, ID.type], as.character(gene.IDs))
    
    # Handle many-to-one warnings
    handle_warning(matched_2)
    
    # Message about the results
    message("Done! ", sum(!is.na(matched)), " genes of ",
            length(matched), " successfully converted.")
    
    # Return results in matrix or vector format
    if (return.Matrix) {
      return(data.frame(GENE.ID = gene.IDs,
                        ENTREZ.ID = GeneNames[matched, "ENTREZID"],
                        stringsAsFactors = FALSE))
    } else {
      return(GeneNames[matched[!is.na(matched)], "ENTREZID"])
    }
    
  } else if (mode == "m2m") {
    # Mouse-to-mouse translation
    if (all(grepl("^ENSM", gene.IDs, perl = TRUE))) {
      ID.type <- "ENSEMBL"
      suppressMessages(GeneNames <- AnnotationDbi::select(
        org.Mm.eg.db::org.Mm.eg.db, gene.IDs, c("ENTREZID", "ENSEMBL"), 
        keytype = "ENSEMBL"))
    } else {
      ID.type <- "SYMBOL"
      suppressMessages(GeneNames <- AnnotationDbi::select(
        org.Mm.eg.db::org.Mm.eg.db, gene.IDs, c("SYMBOL", "ENTREZID"), 
        keytype = "ALIAS"))
    }
    
    # Match results
    matched <- match(as.character(gene.IDs), GeneNames[, ID.type])
    matched_2 <- match(GeneNames[, ID.type], as.character(gene.IDs))
    
    # Handle many-to-one warnings
    handle_warning(matched_2)
    
    # Message about the results
    message("Done! ", sum(!is.na(matched)), " genes of ",
            length(matched), " successfully converted.")
    
    # Return results in matrix or vector format
    if (return.Matrix) {
      return(data.frame(GENE.ID = gene.IDs,
                        ENTREZ.ID = GeneNames[matched, "ENTREZID"],
                        stringsAsFactors = FALSE))
    } else {
      return(GeneNames[matched[!is.na(matched)], "ENTREZID"])
    }
    
  } else if (mode == "m2h") {
    
    # Obtain orthologs
    all.IDs <- babelgene::orthologs(gene.IDs, species = 'mouse', human = FALSE)
    
    # Message about the results
    message("Done! ", length(na.omit(all.IDs$entrez)), 
                      " mouse genes successfully converted to human genes.")
    
    # Return results in matrix or vector format
    if (return.Matrix) {
      return(data.frame(GENE.ID = all.IDs$human_symbol,
                        ENTREZ.ID = all.IDs$human_entrez,
                        stringsAsFactors = FALSE))
    } else {
      return(all.IDs$human_entrez)
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
  
  # Ensure MetaData is loaded
  if (!exists("MetaData", envir = globalenv())) {
    MetaData <- NULL
    data("MetaData", package = "TFEA.ChIP", envir = environment())
  }
  
  # Select the relevant columns from MetaData
  Index <- dplyr::select(MetaData, Accession, TF)
  
  # Apply TF filter if provided
  if (!is.null(TFfilter)) {
    Index <- Index[Index$TF %in% TFfilter, ]
  }
  
  # Apply ENCODE filter if requested
  if (encodeFilter == TRUE) {
    Index <- Index[grepl("^wg|^ENC", Index$Accession), ]
  }
  
  # Check if any results are left, and return or stop with an error
  if (nrow(Index) == 0) {
    stop("No ChIP-Seq dataset in the database follows your conditions.")
  } else {
    return(Index)
  }
}


#### Run TFEA ####

contingency_matrix <- function(test_list, control_list = NULL,
                               chip_index = get_chip_index()) {
  
  #' @title Compute 2x2 Contingency Matrices for ChIP-Seq Enrichment
  #' @description
  #' This function computes 2x2 contingency matrices to assess the overlap between 
  #' a test set of genes and TF binding targets derived from a ChIP-Seq database. 
  #' The matrices are constructed for each ChIP experiment listed in the 
  #' provided index, comparing the number of test and control genes 
  #' that are bound (or not bound) by each transcription factor.
  #' @param test_list A vector of Entrez gene IDs representing the test group 
  #' (e.g., differentially expressed genes).
  #' @param control_list A vector of Entrez gene IDs to be used as the control group. 
  #'   If NULL (default), all genes in the ChIP database not present in `test_list` are used.
  #' @param chip_index A data frame containing metadata for ChIP experiments. Must contain 
  #'   an `Accession` column matching entries in the ChIP database. Defaults to 
  #' the result of `get_chip_index()`.
  #'   `test_list` and `control_list` to simulate a null distribution.
  #' @return A named list of 2x2 contingency matrices, one per ChIP experiment. 
  #'   Each matrix has:
  #'   - Rows: Test group, Control group
  #'   - Columns: Number of genes bound (Positive), and not bound (Negative) by the TF
  #' @examples
  #' data('Genes.Upreg', package = 'TFEA.ChIP')
  #' cm_list <- contingency_matrix(Genes.Upreg)
  #' @export

  
  # Load ChIPDB if not already available
  if (!exists("ChIPDB")) {
    data("ChIPDB", package = "TFEA.ChIP", envir = environment())
  }
  
  # Check data in ChIPDB
  if (!("Gene Keys" %in% names(ChIPDB)) || !("ChIP Targets" %in% names(ChIPDB))) {
    stop("ChIPDB must contain 'Gene Keys' and 'ChIP Targets'.", call. = FALSE)
  }
  
  # Convert ChIPDB to list if it is a matrix
  if (is.matrix(ChIPDB)) {
    ChIPDB <- matrixDB_to_listDB(ChIPDB)
  }
  
  # Generate control gene list if not provided
  if (is.null(control_list)) {
    # Exclude test_list from control
    control_list <- setdiff(ChIPDB[["Gene Keys"]], test_list)  
  } else {
    # Exclude test_list from provided control_list
    control_list <- setdiff(control_list, test_list)  
  }

  # Ensure chip_index is valid
  if (is.null(chip_index) || nrow(chip_index) == 0) {
    stop("The chip_index is empty or does not contain valid accession IDs. Please check the input.")
  }
  
  # Get valid chip accession IDs
  chip_accessions <- intersect(chip_index$Accession, names(ChIPDB[["ChIP Targets"]]))
  
  # Function to compute contingency matrix for each chip
  compute_contingency_matrix <- function(acc, test_genes, control_genes, ChIPDB) {
    chip_targets <- ChIPDB[["Gene Keys"]][ChIPDB[["ChIP Targets"]][[acc]]]
    
    # Calculate the counts for the contingency matrix
    pos_test <- sum(test_genes %in% chip_targets)
    pos_control <- sum(control_genes %in% chip_targets)
    neg_test <- length(test_genes) - pos_test
    neg_control <- length(control_genes) - pos_control
    
    # Create and name the contingency matrix
    cont_matrix <- matrix(c(pos_test, pos_control, neg_test, neg_control),
                          nrow = 2, byrow = TRUE)
    rownames(cont_matrix) <- c("Test", "Control")
    colnames(cont_matrix) <- c("Positive", "Negative")
    
    return(cont_matrix)
  }
  
  # Compute contingency matrices for each chip
  CM_list <- lapply(chip_accessions, compute_contingency_matrix, 
                    test_genes = test_list, control_genes = control_list, 
                    ChIPDB = ChIPDB)
  
  names(CM_list) <- chip_accessions
  
  return(CM_list)
}


getCMstats <- function(CM_list, chip_index = get_chip_index()) {
  
  #' @title Generate statistical parameters from contingency matrices
  #' @description This function computes Fisher's exact test for each matrix
  #' in a list of contingency matrices (e.g., output from `contingency_matrix`).
  #' It returns a data frame containing ChIP-Seq experiment accession IDs,
  #' tested transcription factors, p-values, odds ratios, and adjusted values.
  #' @param CM_list A list of contingency matrices, typically from `contingency_matrix`.
  #' @param chip_index A data frame with ChIP accession IDs and TFs from 
  #' the `get_chip_index` function.
  #' If not provided, the entire internal database is used.
  #' @return A data frame with the ChIP experiment ID, TF, p-value, odds-ratio, 
  #' and other derived statistics.
  #' @export
  #' @examples
  #' data('Genes.Upreg', package = 'TFEA.ChIP')
  #' CM_list_UP <- contingency_matrix(Genes.Upreg)
  #' stats_mat_UP <- getCMstats(CM_list_UP)
  
  # Perform Fisher's exact test on each matrix and calculate OR.SE
  fisher_results <- lapply(seq_along(CM_list), function(i) {
    res <- stats::fisher.test(CM_list[[i]])
    OR <- res$estimate
    OR.SE <- OR * sqrt(1 / CM_list[[i]][1, 1] + 1 / CM_list[[i]][1, 2] + 
                         1 / CM_list[[i]][2, 1] + 1 / CM_list[[i]][2, 2])
    data.frame(OR = OR, p.value = res$p.value, OR.SE = OR.SE)
  })
  fisher_df <- do.call(rbind, fisher_results)
  rownames(fisher_df) <- names(CM_list)
  
  # Ensure chip_index is matched by Accession
  chip_index <- chip_index[match(names(CM_list), chip_index$Accession), ]
  
  # Create statistics data frame
  statMat <- data.frame(
    Accession = chip_index$Accession, 
    TF = chip_index$TF, 
    p.value = fisher_df$p.value, 
    OR = fisher_df$OR, 
    OR.SE = fisher_df$OR.SE, 
    stringsAsFactors = FALSE
  )
  
  # Log transformation of OR and adjustment of p-values
  statMat$log2.OR <- log2(statMat$OR)
  statMat$log2.OR[is.infinite(statMat$log2.OR)] <- NA
  statMat$adj.p.value <- stats::p.adjust(statMat$p.value, "fdr")
  statMat$log10.adj.pVal <- -log10(statMat$adj.p.value)
  statMat$log10.adj.pVal[is.infinite(statMat$log10.adj.pVal)] <- NA
  
  # Cap extreme values for distance calculation
  tmpOR <- statMat$OR
  tmpOR[is.infinite(tmpOR)] <- ifelse(statMat$OR == Inf, 
                                   max(statMat$OR, na.rm = TRUE), 
                                     min(statMat$OR, na.rm = TRUE))
  
  # Calculate Euclidean distance
  statMat$distance <- sapply(seq_along(statMat$Accession), function(i) {
    if (statMat$OR[i] > 1) {
      d <- sqrt(sum((c(statMat$log10.adj.pVal[i], tmpOR[i]) - c(0, 1))^2))
    } else if (statMat$OR[i] <= 1 && statMat$OR[i] > 0) {
      d <- -sqrt(sum((c(statMat$log10.adj.pVal[i], 1 / tmpOR[i]) - c(0, 1))^2))
    } else {
      d <- 0
    }
    return(d)
  })
  
  # Merge with MetaData
  if (!exists("MetaData")) {
    MetaData <- NULL
    data("MetaData", package = "TFEA.ChIP", envir = environment())
  }
  statMat <- merge(MetaData[, c("Accession", "Cell", "Treatment")], 
                       statMat, by = "Accession")
  
  # Sort by Euclidean distance in descending order
  statMat <- statMat[order(statMat$distance, decreasing = TRUE, 
                                                        na.last = TRUE), ]
  
  return(statMat)
}

rankTFs <- function(resultsTable,
                    rankMethod = "gsea", makePlot = FALSE,
                    plotTitle = "TF ranking") {
  
  #' @title Rank the TFs in the output from 'getCMstats'
  #' @description Rank the TFs in the output from 'getCMstats' using
  #' Wilcoxon rank-sum test or a GSEA-like approach.
  #' @param resultsTable Output from the function 'getCMstats'
  #' @param rankMethod "wilcoxon" or "gsea".
  #' @param makePlot (Optional) For rankMethod="gsea". If TRUE, generates 
  #' a plot for TFs with a p-value < 0.05.
  #' @param plotTitle (Optional) Title for the plot.
  #' @return data frame containing:
  #' For Wilcoxon rank-sum test: rank, TF name, test statistic
  #' ('wilc_W), p-value, Freeman's theta, epsilon-squared, and effect size
  #' For GSEA-like ranking: TF name, enrichment score, argument,
  #'   p-value, number of ChIPs
  #' @export rankTFs
  #' @examples
  #' data('Genes.Upreg', package = 'TFEA.ChIP')
  #' CM_list_UP <- contingency_matrix(Genes.Upreg)
  #' stats_mat_UP <- getCMstats(CM_list_UP)
  #' rankTFs(stats_mat_UP)
  
  #### Input format
  rankMethod <- tolower(rankMethod)
  stopifnot(rankMethod %in% c("wilcoxon", "gsea"))
  
  # Check the input type (TFEA.ChIP ORA analysis)
  cols_needed <- c("Accession", "TF", "distance")
  if (any(!cols_needed %in% colnames(resultsTable))) {
    stop("Input error: resultsTable doesn't contain all",
         "the columns required ('Accession', 'TF', 'distance')")
  }
  
  #### GSEA ranking
  if (rankMethod == "gsea") {
    chipSets <- lapply(unique(resultsTable$TF), function(tf) {
      resultsTable$Accession[resultsTable$TF == tf]
    })
    names(chipSets) <- unique(resultsTable$TF)
    
    TFrank <- lapply(chipSets, function(acc_list, statMat) {
      tf <- statMat$TF[statMat$Accession == acc_list[1]]
      res <- GSEA_EnrichmentScore(statMat$Accession, acc_list, 1, 
                                  statMat$distance)
      
      shuffled <- rep( list( statMat$Accession ), 100)
      shuffled <- lapply( shuffled, sample )
      
      shuffled.ES <- sapply(shuffled, function(j, acc, chipDist) {
        tmp.ES <- GSEA_EnrichmentScore(j, acc, 1, chipDist)$ES
        return(tmp.ES)
      }, acc = acc_list, chipDist = statMat$distance)
      
      pVal <- mean(abs(shuffled.ES) > abs(res$ES))
      
      return(data.frame(
        "TF" = tf,
        "ES" = res$ES,
        "arg.ES" = res$arg.ES,
        "pVal" = pVal,
        "numberOfChIPs" = length(acc_list),
        stringsAsFactors = FALSE
      ))
    }, statMat = resultsTable)
    
    TFrank <- do.call(rbind, TFrank)
    
    if (makePlot) {
      plot_df <- TFrank[TFrank$pVal < 0.05 & !is.na(TFrank$pVal), ]
      sub1 <- subset(plot_df, ES > 0)
      sub2 <- subset(plot_df, ES < 0)
      
      if (any(resultsTable$distance == 0)) {
        mid <- max(which(resultsTable$distance == 0)) -
          (min(which(resultsTable$distance == 0)) / 2)
      } else {
        mid <- sum(resultsTable$distance > 0) + 0.5
      }
      
      p <- ggplot2::ggplot() +
        ggplot2::geom_point(
          ggplot2::aes(x = arg.ES, y = ES, color = TF),
          data = plot_df, size = 3, alpha = .5
        ) +
        ggplot2::ylim(-1.5, 1.5) +
        ggplot2::theme_bw() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_point(ggplot2::aes(x = mid, y = 0), color = "black") +  
        ggplot2::theme(legend.position = 'none',
                       plot.title = ggplot2::element_text(hjust = .5, 
                                                          face = "bold", size = 16)) +
        ggplot2::labs(title = plotTitle, 
                      y = "Enrichment Score", 
                      x = "ChIP-seq ranking")
      
      p <- plotly::ggplotly(p)            
      
      
      return(list(TF_ranking = TFrank, TFranking_plot = p))
    } else {
      return(TFrank)
    }
    
    #### Wilcoxon rank-sum test
  } else if (rankMethod == "wilcoxon") {
    chip_ranking <- data.frame(
      acc = resultsTable$Accession,
      TF = resultsTable$TF,
      rank = seq_along(resultsTable$Accession),
      stringsAsFactors = FALSE
    )
    
    TF_wilcox <- lapply(unique(chip_ranking$TF), function(tf, chip_ranking) {
      tmp <- chip_ranking
      tmp$TF[tmp$TF != tf] <- "other"
      tmp$TF <- as.factor(tmp$TF)
      
      wilTest <- stats::wilcox.test(rank ~ TF, data = tmp)
      
      return(data.frame(
        TF = tf,
        wilc_W = wilTest[["statistic"]],
        wilc_pval = wilTest[["p.value"]],
        freeman_theta = rcompanion::freemanTheta(x = tmp$rank, g = tmp$TF),
        epsilon_squared = rcompanion::epsilonSquared(x = tmp$rank, g = tmp$TF),
        wilc_r = rcompanion::wilcoxonR(x = tmp$rank, g = tmp$TF),
        stringsAsFactors = FALSE
      ))
    }, chip_ranking = chip_ranking)
    
    TF_wilcox <- do.call(rbind, TF_wilcox)
    TF_wilcox$rank <- seq_len(nrow(TF_wilcox))
    rownames(TF_wilcox) <- TF_wilcox$rank
    TF_wilcox <- TF_wilcox[, c("rank", "TF", "wilc_W", "wilc_pval",
                               "freeman_theta", "epsilon_squared", "wilc_r")]
    
    return(TF_wilcox)
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
  #' @param correl.vector A vector with the correlations (such as signal to
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
  
  # Ensure correl.vector is provided when weighted score is not 0
  if (weighted.score.type != 0 && is.null(correl.vector)) {
    stop("correl.vector must be provided when weighted.score.type is not 0")
  }
  
  tag.indicator <- as.integer(sign(match(gene.list, gene.set, nomatch = 0)))
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- sum(tag.indicator)  # Number of hits
  Nm <- N - Nh  # Number of misses
  
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  } else {
    if (is.null(correl.vector)) {
      stop("correl.vector must be provided for weighted scores.")
    }
  }
  
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector^alpha)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
  
  if (sum.correl.tag > 0) {
    norm.tag <- 1 / sum.correl.tag
    norm.no.tag <- 1 / Nm
    RES <- cumsum(tag.indicator * correl.vector * norm.tag -
                    no.tag.indicator * norm.no.tag)
    max.ES <- max(RES, na.rm = TRUE)
    min.ES <- min(RES, na.rm = TRUE)
    
    if (abs(max.ES) > abs(min.ES)) {
      ES <- signif(max.ES, digits = 5)
      arg.ES <- match(max.ES, RES)
    } else {
      ES <- signif(min.ES, digits = 5)
      arg.ES <- match(min.ES, RES)
    }
    
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, 
                                indicator = tag.indicator))
  } else {
    # Return NaN or NA values explicitly
    RES <- rep(NaN, length(gene.list))
    indicator <- rep(0, length(gene.list))
    return(list(ES = NA_real_, arg.ES = NA_integer_, 
                                RES = RES, indicator = tag.indicator))
  }
}

GSEA_ESpermutations <- function(gene.list, gene.set, weighted.score.type = 0,
                                correl.vector = NULL, perms = 1000) {
  
  #' @title Calculate enrichment scores for a permutation test.
  #' @description Function to calculate enrichment scores over a randomly 
  #' ordered gene list.
  #' @param gene.list Vector of gene Entrez IDs.
  #' @param gene.set A gene set, e.g. gene IDs corresponding to a ChIP-Seq
  #' experiment's peaks.
  #' @param weighted.score.type Type of score: weight:
  #' 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
  #' @param correl.vector A vector with the correlations (such as signal to
  #' noise scores) corresponding to the genes in the gene list
  #' @param perms Number of permutations
  #' @return Vector of Enrichment Scores for a permutation test.
  #' @examples
  #' GSEA_ESpermutations(gene.list=c('3091','2034','405','55818'),
  #' gene.set=c('2034','112399','405'), perms=10)
  
  # Ensure correl.vector is provided when weighted score is not 0
  if (weighted.score.type != 0 && is.null(correl.vector)) {
    stop("correl.vector must be provided when weighted.score.type is not 0")
  }
  
  tag.indicator <- as.integer(sign(match(gene.list, gene.set, nomatch = 0)))
  N <- length(gene.list)
  Nh <- sum(tag.indicator)  # Number of hits in the gene set
  Nm <- N - Nh  # Number of misses
  
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  } else {
    if (is.null(correl.vector)) {
      stop("correl.vector must be provided for weighted scores.")
    }
  }
  
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector^alpha)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
  
  if (sum.correl.tag > 0) {
    norm.tag <- 1 / sum.correl.tag
    norm.no.tag <- 1 / Nm
    
    pES <- sapply(
      seq_len(perms),
      function(i) {
        mInd <- sample(tag.indicator)
        sum.correl.tag <- sum(correl.vector[mInd == 1])
        mNorm <- 1 / sum.correl.tag
        
        RES <- cumsum(mInd * correl.vector * mNorm - (1 - mInd) * norm.no.tag)
        return(max(abs(RES), na.rm = TRUE))
      }
    )
    return(pES)
    
  } else {
    return(rep(NA_real_, perms))  # Return NA with explicit type
  }
}

GSEA_run <- function(gene.list, LFC, chip_index = get_chip_index(),
                     get.RES = FALSE, RES.filter = NULL, perms = 1000) {
  
  #' @title Function to run a GSEA analysis
  #' @description Analyzes the distribution of TFBS across a sorted list 
  #' of genes.
  #' @param gene.list List of Entrez IDs ordered by their fold change.
  #' @param LFC Vector of log2(Fold Change) values.
  #' @param chip_index Data frame containing accession IDs of ChIPs and the 
  #' tested TFs.
  #' If not provided, the entire internal database will be used.
  #' @param get.RES (Optional) boolean. If TRUE, stores Running Enrichment 
  #' Scores for selected TFs.
  #' @param RES.filter (Optional) character vector. When get.RES == TRUE, 
  #' specifies which TF's RES to store.
  #' @param perms (Optional) integer. Number of permutations for the 
  #' enrichment test.
  #' @return A list containing:
  #' - Enrichment.table: Data frame with enrichment scores and p-values for 
  #' each ChIP-Seq experiment.
  #' - RES (optional): List of running sums for each ChIP-Seq.
  #' - indicators (optional): List of binary vectors indicating matches between
  #' the gene list and gene sets.
  #' @export GSEA_run
  #' @examples
  #' data('hypoxia', package = 'TFEA.ChIP')
  #' hypoxia <- preprocessInputData(hypoxia)
  #' chip_index <- get_chip_index(TFfilter = c('HIF1A', 'EPAS1', 'ARNT'))
  #' GSEA.result <- GSEA_run(hypoxia$Genes, hypoxia$log2FoldChange, chip_index,
  #' get.RES = TRUE)
  
  # Load necessary data if not already loaded
  if (!exists("ChIPDB")) {
    data("ChIPDB", package = "TFEA.ChIP", envir = environment())
  }
  if (is.matrix(ChIPDB)) {
    ChIPDB <- matrixDB_to_listDB(ChIPDB)
  }
  if (!exists("MetaData")) {
    data("MetaData", package = "TFEA.ChIP", envir = environment())
  }
  
  # Initialize result vectors
  enrichmentScore <- numeric()
  pval <- numeric()
  enrichmentArg <- numeric()
  
  res <- if (get.RES) list() else NULL
  ind <- if (get.RES) list() else NULL
  
  # Create progress bar
  pbar <- txtProgressBar(min = 0, max = nrow(chip_index), style = 3)
  
  # Filtered indices for chip_index
  valid_indices <- logical(nrow(chip_index))
  
  for (i in seq_along(chip_index$Accession)) {
    acc <- chip_index$Accession[i]
    targets <- ChIPDB[["Gene Keys"]][ChIPDB[["ChIP Targets"]][[acc]]]
    
    if (sum(targets %in% gene.list) > 10) {
      result <- GSEA_EnrichmentScore(gene.list, targets, 1, LFC)
      shuffled.ES <- GSEA_ESpermutations(gene.list, targets, 1, LFC, perms)
      shuffled.ES <- shuffled.ES[!is.na(shuffled.ES)]
      
      enrichmentScore <- c(enrichmentScore, result$ES)
      pval <- c(pval, sum(shuffled.ES >= abs(result$ES)) / length(shuffled.ES))
      enrichmentArg <- c(enrichmentArg, result$arg.ES)
      
      if (get.RES) {
        if (is.null(RES.filter) || acc %in% RES.filter) {
          res[[acc]] <- result$RES
          ind[[acc]] <- result$indicator
        }
      }
      valid_indices[i] <- TRUE  # Mark valid index
    } else {
      valid_indices[i] <- FALSE  # Mark invalid index
    }
    
    # Update progress bar
    setTxtProgressBar(pbar, i)
  }
  close(pbar)
  
  # Adjust p-values
  pval.adj <- stats::p.adjust(pval, "fdr")
  
  # Create enrichment table
  enrichmentTable <- data.frame(
    Accession = chip_index$Accession[valid_indices],
    TF = chip_index$TF[valid_indices],
    ES = enrichmentScore,
    p.val = pval,
    pval.adj = pval.adj,
    Arg.ES = enrichmentArg,
    stringsAsFactors = FALSE
  )
  
  enrichmentTable <- merge(MetaData[, c("Accession", "Cell", "Treatment")],
                           enrichmentTable, by = "Accession")
  
  # Filter out rows with NA adjusted p-values
  enrichmentTable <- enrichmentTable[!is.na(enrichmentTable$pval.adj), ]
  
  if (get.RES) {
    return(list(Enrichment.table = enrichmentTable, RES = res, indicators = ind))
  } else {
    return(enrichmentTable)
  }
}


#### Plotting functions ####

plot_CM <- function(CM.statMatrix, plot_title = NULL,
                    specialTF = NULL, TF_colors = NULL) {
  
  #' @title Interactive HTML Plot for Transcription Factor Enrichment
  #' @description Generates an interactive HTML plot from a transcription 
  #' factor enrichment table, output of the function 'getCMstats'.
  #' @param CM.statMatrix Output of the function 'getCMstats', a data frame 
  #' containing Accession ID, Transcription Factor, Odds Ratio, p-value, and 
  #' adjusted p-value.
  #' @param plot_title The title for the plot (default: "Transcription Factor 
  #' Enrichment").
  #' @param specialTF (Optional) Named vector of TF symbols to be highlighted 
  #' in the plot, allowing for grouped color representation.
  #' @param TF_colors (Optional) Colors to highlight TFs specified in 
  #' specialTF.
  #' @return A plotly scatter plot.
  #' @export plot_CM
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("The 'plotly' package is required for this function to work. Please install it.",
                                                                            call. = FALSE)
  }
  
  # Set default plot title
  plot_title <- ifelse(is.null(plot_title), 
                        "Transcription Factor Enrichment", plot_title)
  
  # Prepare highlighting for transcription factors
  if (is.null(specialTF)) {
    CM.statMatrix$highlight <- "Others"
    TF_colors <- "grey"
  } else {
    if (is.null(TF_colors)) {
      TF_colors <- c(setNames(RColorBrewer::brewer.pal(min(length(specialTF),
                                                        8), "Set1"), specialTF), "grey")
    } else {
      # Ensure colors are provided for all special TFs, default to grey for others
      TF_colors <- c(TF_colors, "grey")
    }
    CM.statMatrix$highlight <- ifelse(CM.statMatrix$TF %in%
                                                specialTF, CM.statMatrix$TF, "Others")
  }
  
  # Load metadata if it does not exist in the environment
  if (!exists("MetaData", envir = .GlobalEnv)) {
    data("MetaData", package = "TFEA.ChIP", envir = .GlobalEnv)
  }
  
  # Merge CM.statMatrix with MetaData
  CM.statMatrix <- merge(CM.statMatrix, MetaData[, c("Accession", "Treatment", "Cell")], 
                                                           by = "Accession", all.x = TRUE)
  colnames(CM.statMatrix)[grepl('Cell.x', colnames(CM.statMatrix))] <- 'Cell'
  colnames(CM.statMatrix)[grepl('Treatment.x', colnames(CM.statMatrix))] <- 'Treatment'
  
  # Replace Inf and -Inf values in Odds Ratio
  CM.statMatrix$OR[CM.statMatrix$OR == Inf] <- max(CM.statMatrix$OR[CM.statMatrix$OR != Inf], 
                                                                                  na.rm = TRUE)
  CM.statMatrix$OR[CM.statMatrix$OR == -Inf] <- min(CM.statMatrix$OR[CM.statMatrix$OR != -Inf], 
                                                                                  na.rm = TRUE)
  
  # Handle p-values of zero for log scale plotting
  if (any(CM.statMatrix$adj.p.value == 0)) {
    finMax <- max(CM.statMatrix$log.adj.pVal[CM.statMatrix$p.value != 0], na.rm = TRUE)
    CM.statMatrix$log.adj.pVal[CM.statMatrix$p.value == 0] <- finMax
  }
  
  # Add Cell and Treatment to the hover text
  CM.statMatrix$pointText <- paste0("<b>Accession:</b> ", CM.statMatrix$Accession, 
                                    "<br><b>TF:</b> ", CM.statMatrix$TF, 
                                    "<br><b>Cell:</b> ", CM.statMatrix$Cell, 
                                    "<br><b>Treatment:</b> ", CM.statMatrix$Treatment, 
                                    "<br><b>Odds Ratio:</b> ", round(CM.statMatrix$OR, 2), 
                                    "<br><b>Adj p-value:</b> ", signif(CM.statMatrix$adj.p.value, 3))
  
  # Plot using ggplot2
  p <- ggplot2::ggplot(CM.statMatrix, ggplot2::aes(x = log10.adj.pVal, 
                                                   y = log2.OR,
                                                   label = ifelse(highlight != "Others", highlight, NA),
                                                   text = pointText)) +
    # Plot grey points first (for "Others" category)
    ggplot2::geom_point(data = CM.statMatrix[CM.statMatrix$highlight == "Others", ], 
                        ggplot2::aes(color = highlight), 
                        size = 2, alpha = 0.7) + 
    # Then plot the colored points on top (for "specialTF")
    ggplot2::geom_point(data = CM.statMatrix[CM.statMatrix$highlight != "Others", ], 
                        ggplot2::aes(color = highlight), 
                        size = 4, alpha = 1, show.legend = TRUE) + 
    ggplot2::scale_color_manual(values = TF_colors) + 
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme_bw() +
    ggplot2::labs(color = "TF") +  # Change legend name to "TF"
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, face = "bold", size = 16))
  
  # Convert to interactive plot using plotly and set tooltips
  p <- plotly::ggplotly(p, tooltip = "text")
  
  return(p)
}


plot_ES <- function(GSEA_result, LFC, plot_title = NULL, specialTF = NULL,
                    Accession = NULL, TF = NULL) {
  
  #' @title Plots Enrichment Score from the output of GSEA.run.
  #' @description Function to plot the Enrichment Score of every member of
  #' the ChIPseq binding database.
  #' @param GSEA_result Returned by GSEA_run.
  #' @param LFC Vector with log2(Fold Change) of every gene with an Entrez ID, 
  #' ordered from highest to lowest.
  #' @param plot_title (Optional) Title for the plot.
  #' @param specialTF (Optional) Named vector of transcription factors (TF)
  #'  to highlight in the plot.
  #' @param Accession (Optional) Vector of dataset IDs to restrict the plot to.
  #' @param TF (Optional) Vector of transcription factor names to restrict 
  #' the plot to.
  #' @return Plotly object combining scatter plot of enrichment scores and a 
  #' log2(fold change) heatmap.
  #' @export
  
  # Extract enrichment table from GSEA result
  enrichTab <- if (is.data.frame(GSEA_result)) {
    GSEA_result
  } else if ("Enrichment.table" %in% names(GSEA_result)) {
    GSEA_result[["Enrichment.table"]]
  } else {
    stop("Invalid GSEA_result format: missing Enrichment.table.")
  }
  
  # Filter based on provided Accession or TF
  if (!is.null(Accession)) {
    enrichTab <- enrichTab[enrichTab$Accession %in% Accession, ]
  }
  if (!is.null(TF)) {
    enrichTab <- enrichTab[enrichTab$TF %in% TF, ]
  }
  
  # Set plot title if not provided
  if (is.null(plot_title)) {
    plot_title <- "Transcription Factor Enrichment"
  }
  
  # Handle special TFs for highlighting
  if (is.null(specialTF)) {
    enrichTab$highlight <- "Other"
    markerColors <- c("Other" = "azure4")
  } else {
    # Function to mark and color the special TFs
    color_list <- highlight_TF(enrichTab, 4, specialTF)
    enrichTab$highlight <- color_list[[1]]
    markerColors <- color_list[[2]]
  }
  
  # Add significance labels
  enrichTab$symbol <- ifelse(enrichTab$pval.adj <= 0.05, 
                             "pVal<=0.05", "pVal>0.05")
  
  # Load metadata for treatment and cell information if available
  if (!exists("MetaData", where = environment())) {
    data("MetaData", package = "TFEA.ChIP", envir = environment())
  }
  
  # Prepare data for ggplot
  enrichTab$Treatment <- MetaData$Treatment[match(enrichTab$Accession, MetaData$Accession)]
  enrichTab$Cell <- MetaData$Cell[match(enrichTab$Accession, MetaData$Accession)]
  
  # Create scatter plot using ggplot2
  es_plot <- ggplot2::ggplot(enrichTab, ggplot2::aes(x = Arg.ES, y = ES, 
                                                     color = highlight, 
                                                     shape = symbol,
                                                     text = paste0(
                                                       Accession,
                                                       "<br>ES:", round(ES, 3),
                                                       "<br>Adj.P.value:", round(pval.adj, 3),
                                                       "<br>Treatment:", Treatment,
                                                       "<br>Cell:", Cell))) +
    ggplot2::geom_point(size = 3) +  
    ggplot2::scale_color_manual(values = markerColors) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = plot_title, x = "Rank", y = "Enrichment Score (ES)") +
    ggplot2::scale_shape_manual(values = c("pVal<=0.05" = 17, "pVal>0.05" = 19)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Prepare log2(Fold Change) data for the bar plot
  lfc_df <- data.frame(x = seq_along(LFC), y = LFC)
  lfc_bar <- ggplot2::ggplot(lfc_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bar(stat = "identity", color = "grey40") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Rank", y = "Log2(FC)") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Convert ggplot objects to plotly for interactive visualization
  es_plotly <- plotly::ggplotly(es_plot, tooltip = c("text"))
  
  lfc_plotly <- plotly::ggplotly(lfc_bar)
  
  # Combine both plots using plotly::subplot
  combined_plot <- plotly::subplot(es_plotly, lfc_plotly, nrows = 2, shareX = TRUE,
                                   heights = c(0.8, 0.2), titleY = TRUE)
  
  return(combined_plot)
}

highlight_TF <- function(enrichTab, column, specialTF) {
  
  #' @title Highlight special TFs in enrichment table.
  #' @description Function to assign special TFs to specific colors and 
  #' highlight them in the enrichment table.
  #' @param enrichTab Data frame containing the enrichment results.
  #' @param column The column index or name of the enrichment table used 
  #' for matching TFs.
  #' @param specialTF A named vector of transcription factors to highlight
  #'  in the plot.
  #' @return A list containing the updated highlight column and the color 
  #' mapping for each TF.
  #' @export
  
  enrichTab$highlight <- "Other"
  
  # Ensure valid matching and assignment
  match_idx <- which(enrichTab[[column]] %in% specialTF)
  if (length(match_idx) > 0) {
    enrichTab$highlight[match_idx] <- enrichTab[[column]][match_idx]
  }
  
  # Generate colors for each special TF
  colors.plot <- rainbow(length(unique(specialTF)))
  
  # Map the TFs to colors
  markerColors <- c(setNames(colors.plot, specialTF), "Other" = "azure4")
  
  return(list(enrichTab$highlight, markerColors))
}


plot_RES <- function(GSEA_result, LFC, plot_title = NULL, line.colors = NULL,
                     line.styles = NULL, Accession = NULL, TF = NULL) {
  
  #' @title Plot Running Enrichment Scores (RES) and Log2 Fold Change (LFC)
  #' @description This function plots the running enrichment scores (RES) from 
  #' a GSEA result, with an additional bar plot showing the log2 fold change 
  #' (LFC) of genes. The RES plot can be filtered by TFs or accession IDs.
  #' @param GSEA_result List returned by GSEA_run, containing the enrichment 
  #' table and RES.
  #' @param LFC Numeric vector containing the log2(Fold Change) of every gene 
  #' with an Entrez ID, ordered from highest to lowest.
  #' @param plot_title (Optional) String specifying the title for the plot.
  #'  Default is "Transcription Factor Enrichment".
  #' @param line.colors (Optional) Character vector specifying colors for each 
  #' line in the RES plot. If NULL, default colors will be used.
  #' @param line.styles (Optional) Character vector specifying line styles for 
  #' each RES line. Possible values are 'solid', 'dash', or 'longdash'.
  #' @param Accession (Optional) Character vector specifying accession IDs to
  #' restrict the plot to. If NULL, all accession IDs will be plotted.
  #' @param TF (Optional) Character vector specifying TF names to restrict 
  #' the plot to. If NULL, all transcription factors will be plotted.
  #' @return A Plotly object containing two subplots: the top one showing the 
  #' running enrichment scores (RES) for the filtered accession IDs or 
  #' TFs, and the bottom one displaying the log2 fold
  #' change (LFC) as a bar plot.
  #' @export plot_RES
  #' @examples
  #' data('GSEA_result', 'log2.FC', package = 'TFEA.ChIP')
  #' plot_RES(GSEA_result, log2.FC, 
  #'          TF = c('E2F4', 'E2F1'),
  #'          Accession = c('ENCSR000DYY.E2F4.GM12878', 
  #'                        'ENCSR000EVJ.E2F1.HeLa-S3'))
  
  # Load MetaData if not already available
  if (!exists("MetaData")) {
    data("MetaData", package = "TFEA.ChIP", envir = environment())
  }
  
  # Extract enrichment table and running sums from GSEA_result
  enrichTab <- GSEA_result[["Enrichment.table"]]
  rSums <- GSEA_result[["RES"]]
  
  # Apply Accession or TF filters
  if (!is.null(Accession) | !is.null(TF)) {
    if (is.null(Accession)) {
      Accession <- enrichTab$Accession[enrichTab$TF %in% TF]
    } 
    if (is.null(TF)) {
      TF <- enrichTab$TF[enrichTab$Accession %in% Accession]
    }
    enrichTab <- enrichTab[match(Accession, enrichTab$Accession), ]
    rSums <- rSums[match(Accession, names(rSums))]
  } else {
    Accession <- enrichTab$Accession
  }
  
  # Set default line parameters
  line.colors <- ifelse(is.null(line.colors), 
                        c("red", "blue", "green", "hotpink", "cyan", 
                          "greenyellow", "gold", "darkorchid", 
                          "chocolate1", "black", "lightpink", 
                          "seagreen")[1:length(Accession)], 
                        line.colors)
  
  plot_title <- ifelse(is.null(plot_title), 
                       "Transcription Factor Enrichment", 
                       plot_title)
  
  # Prepare data for ggplot
  MetaData <- MetaData[match(Accession, MetaData$Accession), ]
  
  # Create a long-format data frame
  plotTab <- data.frame(
    x = unlist(lapply(rSums, seq_along)),  
    RES = unlist(rSums),                  # RES values
    Accession = rep(Accession, lengths(rSums)), 
    TF = rep(MetaData$TF, lengths(rSums)), 
    stringsAsFactors = FALSE
  )
  
  # Create ggplot with multiple RES lines
  res_plot <- ggplot2::ggplot(plotTab, 
                  ggplot2::aes(x = x, y = RES, color = Accession)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_hline(yintercept = 0, color = "grey30", size = .5) + 
    ggplot2::theme_bw() +
    ggplot2::labs(title = plot_title, x = "Rank", y = "Enrichment Score (RES)") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16), 
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank()   
    )

  # Plot the LFC as a bar plot underneath
  lfc_df <- data.frame(x = seq_along(LFC), y = LFC)
  lfc_bar <- ggplot2::ggplot(lfc_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bar(stat = "identity",  color = "grey40") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "Log2(FC)") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16),  
      panel.grid.major = ggplot2::element_blank(),  
      panel.grid.minor = ggplot2::element_blank()  
    )
  
  # Convert to plotly
  res_plotly <- plotly::ggplotly(res_plot)
  lfc_plotly <- plotly::ggplotly(lfc_bar)
  
  # Combine the plots 
  combined_plot <- plotly::subplot(res_plotly, lfc_plotly, 
                                   nrows = 2, shareX = TRUE,
                                   heights = c(0.8, 0.2), titleY = TRUE)
  
  return(combined_plot)
}

metaanalysis_fx <- function(dat) {
  #' @title Perform Meta-analysis for each TF
  #' @description Conducts a random-effects meta-analysis of odds ratios (OR) 
  #' and standard errors (OR.SE) for each TF using the `meta` package.
  #' @param dat A data frame with columns: TF, OR, OR.SE, Accession, adj.pval.
  #' @return A list with:
  #'   - summary: a data frame of ranked meta-analysis results per TF
  #'   - results: a named list of raw meta-analysis objects from the 
  #'    `meta` package
  #' @export
  #' @examples
  #' df <- data.frame(TF = c('A', 'A', 'A', 'B', 'B'),
  #'                  OR = c(1, 1.2, 1.23, 4, 4.5),
  #'                  OR.SE = c(1e-5, 5e-4, 2e-4, 1e-3, 1e-2),
  #'                  Accession = c('Chip1', 'Chip2', 'Chip3', 'Chip4', 'Chip5'),
  #'                  adj.pval = c(1e-5, 5e-4, 2e-4, 1e-3, 1e-2))
  #' res <- metaanalysis_fx(df)
  
  error_tfs <- character()
  tf_list <- unique(dat$TF)
  meta_results <- list()
  
  result_df <- purrr::map_dfr(tf_list, function(tf) {
    tf_data <- dplyr::filter(dat, TF == tf)
    
    model <- tryCatch({
      meta::metagen(
        TE = tf_data$OR,
        seTE = tf_data$OR.SE,
        data = tf_data,
        sm = "MD",
        common = FALSE,
        random = TRUE,
        method.tau = "REML",
        cluster = tf_data$Accession,
        studlab = tf_data$Accession,
        title = tf,
        prediction = TRUE
      )
    }, error = function(e) {
      error_id <- limma::strsplit2(tf_data$Accession[1], "[.]")[, 2]
      error_tfs <<- c(error_tfs, error_id)
      return(NULL)
    })
    
    if (!is.null(model)) {
      meta_results[[tf]] <<- model
      
      tibble::tibble(
        TF = tf,
        k = model$k,
        OR = model$TE.random,
        OR.SE = model$seTE.random,
        CI95lower = model$lower.random,
        CI95upper = model$upper.random,
        pval = model$pval.random
      )
    }
  })
  
  if (nrow(result_df) == 0) {
    stop("No valid meta-analysis results were produced.")
  }
  
  result_df <- result_df %>%
    dplyr::mutate(
      pval.adj = stats::p.adjust(pval, method = "fdr"),
      pval.adj.safe = ifelse(pval.adj == 0, .Machine$double.xmin, pval.adj),
      rank_score = abs(log2(OR)) * -log10(pval.adj.safe)
    ) %>%
    dplyr::arrange(dplyr::desc(rank_score)) %>%
    dplyr::select(-pval.adj.safe)
  
  if (length(error_tfs) > 0) {
    message("Meta-analysis failed for the following TFs: ",
            paste(unique(error_tfs), collapse = ", "),
            ". Please check their individual ChIP-seq values.")
  }
  
  return(list(summary = result_df, results = meta_results))
}


filter_expressed_TFs <- function(Table, chip_index, TFfilter = NULL, encodeFilter = FALSE) {
  
  #' @title Filter Expressed TFs
  #' @description Filters TFs based on their expression status in the input dataset.
  #' This function identifies expressed TFs by intersecting the input gene list
  #' with the `chip_metadata` dataset.
  #' @param Table A data frame containing gene expression data with a `Genes` column.
  #' @param chip_index A data frame containing ChIP-Seq dataset accession IDs 
  #' and associated TFs.
  #' @param TFfilter (Optional) A character vector of TFs to filter.
  #' @param encodeFilter (Optional) Logical; if TRUE, applies ENCODE filtering 
  #' to ChIP-Seq data.
  #' @return A filtered `chip_index` data frame containing only expressed TFs.
  #' @export
  
  # Load chip_metadata if not already available
  if (!exists("chip_metadata")) {
    data("chip_metadata", package = "TFEA.ChIP", envir = environment())
  }
  
  # Identify expressed genes
  expressed_genes <- intersect(Table$Genes, chip_metadata$EntrezID)
  
  # Identify removed genes
  rm_genes <- setdiff(chip_metadata$EntrezID, Table$Genes)
  rm_genes_names <- chip_metadata %>%
    filter(EntrezID %in% rm_genes) %>%
    pull(tf.name)
  
  # Provide informative messages about removed genes
  if (length(rm_genes_names) > 0) {
    message("The following TFs were excluded from the analysis because they were not found in the provided dataset and are considered not expressed. To modify this behavior, set `expressed` to FALSE:\n", 
            paste(sort(unique(rm_genes_names)), collapse = ", "))
  } else {
    message("No genes were removed from the analysis as all genes were present in the provided dataset.")
  }
  
  # Filter TFs based on expressed genes
  filtered_tfs <- chip_metadata %>%
    filter(EntrezID %in% expressed_genes) %>%
    pull(chip.name)
  
  # Apply additional TF filtering if provided
  if (!is.null(TFfilter)) {
    filtered_tfs <- intersect(filtered_tfs, TFfilter)
  }
  
  # Update chip index with filtered TFs
  chip_index <- get_chip_index(TFfilter = filtered_tfs, encodeFilter = encodeFilter)
  
  return(chip_index)
}

analysis_from_table <- function(inputData, mode = "h2h", 
                                interest_min_LFC = -Inf, 
                                interest_max_LFC = Inf,
                                control_min_LFC = -0.25, 
                                control_max_LFC = 0.25,
                                interest_min_pval = 0, 
                                interest_max_pval = 0.05,
                                control_min_pval = 0.5, 
                                control_max_pval = 1, 
                                expressed = TRUE,
                                encodeFilter = FALSE,
                                TFfilter = NULL,
                                method = "ora") {
  
  #' @title Analysis from Input Table
  #' @description Performs gene expression analysis, filtering genes and
  #' TFs based on specified thresholds. It calculates statistics using 
  #' overrepresentation analysis (ORA) or gene set enrichment analysis (GSEA).
  #'
  #' @param inputData A data frame containing gene expression data.
  #' @param mode Character string specifying the mode: 'h2h', 'm2h', 'm2m'.
  #' @param interest_min_LFC Minimum LFC for selected genes of interest.
  #' @param interest_max_LFC Maximum LFC for selected genes of interest.
  #' @param control_min_LFC Minimum LFC for control genes.
  #' @param control_max_LFC Maximum LFC for control genes.
  #' @param interest_min_pval Minimum p-value for genes of interest.
  #' @param interest_max_pval Maximum p-value for genes of interest.
  #' @param control_min_pval Minimum p-value for control genes.
  #' @param control_max_pval Maximum p-value for control genes.
  #' @param expressed Logical; filter TFs expressed in input data.
  #' @param encodeFilter Logical; apply ENCODE filtering to ChIP-seq data.
  #' @param TFfilter Character vector of transcription factors to filter (optional).
  #' @param method Analysis method: 'ora' (overrepresentation) or 'gsea' (gene set enrichment).
  #' @return A matrix with calculated statistics (e.g., p-values, odds ratios).
  #' @export
  #' @examples
  #' data('hypoxia_DESeq',package='TFEA.ChIP')
  #' res <- analysis_from_table(hypoxia_DESeq, interest_min_LFC = 1)
  
  # Validate mode
  if (!mode %in% c("h2h", "m2h", "m2m")) {
    stop("Invalid mode. Use 'h2h', 'm2h', or 'm2m'.")
  }
  
  # Preprocess input data
  cat("Preprocessing input data...\n")
  Table <- preprocessInputData(inputData, mode)
  if (is.null(Table)) stop("Error: Failed to preprocess input data.")

  # Load ChIP-seq data if not already loaded
  if (!exists("ChIPDB")) {
    data("ChIPDB", package = "TFEA.ChIP", envir = environment())
  }
  if (!exists("chip_metadata")) {
    data("chip_metadata", package = "TFEA.ChIP", envir = environment())
  }
  
  # Filter TFs
  if (!is.null(TFfilter)) {
    TFfilter <- chip_metadata %>%
      filter(tf.name %in% TFfilter) %>%
      pull(chip.name)
  }
  
  # Retrieve chip index based on filters
  cat("Retrieving ChIP index...\n")
  chip_index <- get_chip_index(TFfilter = TFfilter, encodeFilter = encodeFilter)
  
  # Filter for expressed TFs if required
  if (expressed) {
    cat("Filtering for expressed transcription factors...\n")
    chip_index <- filter_expressed_TFs(Table, chip_index, TFfilter, encodeFilter)
  }
  
  # Perform analysis
  cat("Performing analysis using the selected method...\n")
  result <- NULL
  if (method == "ora") {
    
    # Select genes based on thresholds
    Genes.Selected <- Select_genes(GeneExpression_df = Table, 
                                   min_LFC = interest_min_LFC, 
                                   max_LFC = interest_max_LFC, 
                                   min_pval = interest_min_pval, 
                                   max_pval = interest_max_pval)
    
    Genes.Control <- Select_genes(GeneExpression_df = Table, 
                                  min_LFC = control_min_LFC, 
                                  max_LFC = control_max_LFC, 
                                  min_pval = control_min_pval, 
                                  max_pval = control_max_pval)
    
    # Warn if no genes are selected
    if (length(Genes.Selected) == 0) {
      warning("No genes selected with the specified thresholds.")
    }
    if (length(Genes.Control) == 0) {
      warning("No control genes found with the specified thresholds.")
    }
    
    CM_list <- contingency_matrix(Genes.Selected, Genes.Control, chip_index)
    if (length(CM_list) == 0) stop("Error: No contingency matrices generated.")
    
    cat("Calculating statistics from contingency matrices...\n")
    result <- getCMstats(CM_list, chip_index = chip_index)
    
  } else if (method == "gsea") {
    
    cat("Running Gene Set Enrichment Analysis ...\n")
    gsea_res <- GSEA_run(Table$Genes, Table$log2FoldChange, 
                       chip_index = chip_index, get.RES = TRUE)
    
    result <- list(result = gsea_res, processed_table = Table)
  } else {
    stop("Invalid method. Use 'ora' or 'gsea'.")
  }
  
  cat("Analysis completed!\n")
  return(result = result)
}



