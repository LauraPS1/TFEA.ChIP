# Load the required package
library(RUnit)
library(TFEA.ChIP)

# Testing units for the TFEA.ChIP package

test_GeneID2Entrez<-function(){
    # Checking symbol translation
    RUnit::checkEquals(
        GeneID2entrez(c("HIF1A","EPAS1")),
        c("3091","2034"))
    # Checking lowercase handling
    RUnit::checkEquals(
        GeneID2entrez("HIF1A"),
        GeneID2entrez("hif1a"))
    # Checking Ensembl translation
    RUnit::checkEquals(
        GeneID2entrez(c("ENSG00000100644","ENSG00000116016")),
        c("3091","2034"))
    # Checking matrix output and NA return for invalid names
    RUnit::checkException(
        is.na(GeneID2entrez("potato",return.Matrix = TRUE)$ENTREZ.ID))
}

test_preprocessInputData<-function(){
    data("hypoxia_DESeq","hypoxia",package = "TFEA.ChIP")
    # Checking output size
    RUnit::checkEquals(
        suppressWarnings(nrow(preprocessInputData(hypoxia))),
        suppressWarnings(length(GeneID2entrez(hypoxia$Genes))))

    RUnit::checkEquals(
        suppressWarnings(nrow(preprocessInputData(hypoxia_DESeq))),
        suppressWarnings(length(GeneID2entrez(rownames(hypoxia_DESeq)))))

    # Cheking every gene in output has been properly translated
    RUnit::checkTrue(!(FALSE %in% (
        suppressWarnings(preprocessInputData(hypoxia)$Genes) %in%
        suppressWarnings(GeneID2entrez(hypoxia$Genes)))))
    RUnit::checkTrue(!(FALSE %in% (
        suppressWarnings(preprocessInputData(hypoxia_DESeq)$Genes) %in%
        suppressWarnings(GeneID2entrez(rownames(hypoxia_DESeq))))))

}

test_Select_genes<-function(){
    data("hypoxia",package="TFEA.ChIP")
    # All differentially expressed genes
    RUnit::checkEquals(
        length( Select_genes( hypoxia ) ),
        dim( hypoxia[ hypoxia$pval.adj<=0.05, ] )[1] )
    # Downregulated genes
    RUnit::checkEquals(
        length( Select_genes( hypoxia, max_LFC = 0 ) ),
        dim( hypoxia[ hypoxia$pval.adj<=0.05 & hypoxia$log2FoldChange<0, ] )[1] )
    # Upregulated genes
    RUnit::checkEquals(
        length( Select_genes( hypoxia, min_LFC = 0 ) ),
        dim( hypoxia[ hypoxia$pval.adj<=0.05 & hypoxia$log2FoldChange>0, ] )[1] )
    # Unresponsive genes
    RUnit::checkEquals(
        length( Select_genes( hypoxia,
            min_pval = 0.5, max_pval = 1, min_LFC = -0.25, max_LFC = 0.25)),
        dim( hypoxia[ hypoxia$pval.adj >= 0.5 &
            hypoxia$log2FoldChange > (-0.25) & hypoxia$log2FoldChange < 0.25, ] )[1] )
}

test_get_chip_index<-function(){
    data("MetaData",package = "TFEA.ChIP")

    # Full database
    RUnit::checkEquals(
        nrow(get_chip_index()),
        nrow(MetaData))
    # Datasets from Encode only
    RUnit::checkEquals(
        nrow(get_chip_index(encodeFilter = TRUE)),
        nrow(MetaData[grepl("^ENC",MetaData$Accession),]))
    # TF "JUND" only
    RUnit::checkEquals(
        nrow(get_chip_index(TFfilter = "JUND")),
        nrow(MetaData[MetaData$TF=="JUND",]))
}

test_contingency_matrix<-function(){
    data("Genes.Upreg",package = "TFEA.ChIP")
    data("ChIPDB",package = "TFEA.ChIP")
    # Sum of every contingency matrix. Since no control gene list
    # is provided, all genes not included in the test list will
    # be used as control.
    cont_mat<-contingency_matrix(Genes.Upreg)
    RUnit::checkTrue(!(FALSE %in%
        (sapply(cont_mat,sum) == nrow(ChIPDB)))
    )
    # Size of the output when selecting ChIP-Seq datasets
    chip_index<-get_chip_index(TFfilter = "JUN")
    RUnit::checkEquals(
        length(contingency_matrix(Genes.Upreg,chip_index = chip_index)),
        nrow(chip_index)
    )
    RUnit::checkEquals(
        names(contingency_matrix(Genes.Upreg,chip_index = chip_index)),
        chip_index$Accession
    )

}

test_getCMstats<-function(){
    data('Genes.Upreg',package = 'TFEA.ChIP')
    CM_list <- contingency_matrix(Genes.Upreg)
 
    # Checking output size
    RUnit::checkEquals(
        nrow(getCMstats(CM_list)),
        length(CM_list))
    # Checking output values
    RUnit::checkEquals(
        getCMstats(CM_list[1])$p.value,
        stats::fisher.test(CM_list[[1]])$p.value)
    RUnit::checkEquals(
        getCMstats(CM_list)$adj.p.value,
        p.adjust(getCMstats(CM_list)$p.value,"fdr"))
}

test_GSEA_EnrichmentScore<-function(){

    # Checking matches
    RUnit::checkEquals(
        GSEA_EnrichmentScore(c(1:10),c(1:5))$indicator,
        c(rep(1,5),rep(0,5)))
    # Checking score and argument
    RUnit::checkTrue(GSEA_EnrichmentScore(c(1:10),c(1:5))$ES<=1)
    RUnit::checkTrue(
        GSEA_EnrichmentScore(c(1:10),c(1:5))$arg.ES<=10)
}

test_txt2GR<-function(){
    data("ARNT.peaks.bed","ARNT.metadata",package = "TFEA.ChIP")
    # Checking output size
    RUnit::checkEquals(
        length(ARNT.gr<-txt2GR(ARNT.peaks.bed,"macs1.4",ARNT.metadata)),
        nrow(ARNT.peaks.bed[ARNT.peaks.bed$X.10.log10.pvalue.<50]))
    # Checking that missing input variables produces errors
    RUnit::checkException(
        ARNT.gr<-txt2GR(ARNT.peaks.bed,fileMetaData = ARNT.metadata))
    RUnit::checkException(
        ARNT.gr<-txt2GR(ARNT.peaks.bed,"macs1.4"))
    # Checking incorrect input variables produces errors
    RUnit::checkException(
        ARNT.gr<-txt2GR(ARNT.peaks.bed,"macs1.4",ARNT.metadata[,1:3]))

}

test_makeChIPGeneDB<-function(){
    data("DnaseHS_db","gr.list", package="TFEA.ChIP")
    # Using this toy datasets, the expected result is a
    # list of one gene set with 54 gene IDs
    
    ChIPDB <- makeChIPGeneDB( DnaseHS_db, gr.list )
    
    RUnit::checkEquals( length( ChIPDB[["ChIP Targets"]] ),1)
    RUnit::checkEquals(
        length(ChIPDB[["ChIP Targets"]][[ 1 ]]),
        54 )
}

test_GSEA_run<-function(){
    data("hypoxia",package = "TFEA.ChIP")
    hypoxia<-preprocessInputData(hypoxia)
    # Selecting a small number of ChIP-Seq datasets to save time
    chip_index<-get_chip_index(TFfilter = c("EPAS1","ARNT"))
    gsea_result <- GSEA_run(hypoxia$Genes, hypoxia$log2FoldChange, chip_index)
    # Checking enrichment table
    RUnit::checkTrue(
        all( gsea_result$Accession %in%
        chip_index$Accession))

    # Checking all ChIPs done in hypoxic conditions are enriched
    RUnit::checkTrue(
        all(gsea_result[gsea_result$Treatment!="none", "ES"] > 0),
        all(gsea_result[
            gsea_result$Treatment!="none", "Arg.ES"] < dim(hypoxia)[1]/2)
    )
    # Checking output with RES and indicator
    RUnit::checkTrue(is.list(
        GSEA_run(hypoxia$Genes, hypoxia$log2FoldChange, chip_index[1:2,],get.RES = TRUE)))


}



test_get_LFC_bar <- function() {
    # Test with normal data
    logFC <- c(0.5, -0.2, 1.5, -0.5, 0.0)
    RUnit::checkTrue(is(plotly::plotly_build(get_LFC_bar(logFC)), "plotly"))
    
    # Test with all positive values
    logFC_pos <- c(1, 2, 3)
    RUnit::checkTrue(is(plotly::plotly_build(get_LFC_bar(logFC_pos)), "plotly"))
    
    # Test with all negative values
    logFC_neg <- c(-1, -2, -3)
    RUnit::checkTrue(is(plotly::plotly_build(get_LFC_bar(logFC_neg)), "plotly"))

    # Test with mixed values
    logFC_mixed <- c(-1, 0, 1, 2)
    RUnit::checkTrue(is(plotly::plotly_build(get_LFC_bar(logFC_mixed)), "plotly"))

    # Test with empty input
    RUnit::checkException(get_LFC_bar(c()))

    # Test with invalid input (non-numeric)
    RUnit::checkException(get_LFC_bar(c("a", "b", "c")))
}


test_meta_analysis_fx <- function() {
    # Test with valid data
    data <- data.frame(
        TF = c('TF1', 'TF1', 'TF2', 'TF2'),
        OR = c(1.5, 1.2, 0.9, 1.1),
        OR.SE = c(0.2, 0.3, 0.1, 0.25),
        Accession = c('D1', 'D2', 'D3', 'D4'),
        adj.pval = c(0.01, 0.02, 0.05, 0.03)
    )
    result <- metaanalysis_fx(data)
    RUnit::checkEquals(nrow(result), 2)

    # Test with edge case: no datasets
    empty_data <- data.frame(TF = character(0), OR = numeric(0), OR.SE = numeric(0), Accession = character(0), adj.pval = numeric(0))
    RUnit::checkException(metaanalysis_fx(empty_data))

    # Test with invalid data (non-numeric in OR or OR.SE)
    invalid_data <- data.frame(
        TF = c('TF1', 'TF1'),
        OR = c('a', 'b'),
        OR.SE = c(0.2, 0.3),
        Accession = c('D1', 'D2'),
        adj.pval = c(0.01, 0.02)
    )
    RUnit::checkException(metaanalysis_fx(invalid_data))

    # Test for results with no valid OR
    no_valid_or_data <- data.frame(
        TF = c('TF1', 'TF1', 'TF2', 'TF2'),
        OR = c(1, 1, 1, 1),
        OR.SE = c(0, 0, 0, 0),
        Accession = c('D1', 'D2', 'D3', 'D4'),
        adj.pval = c(1, 1, 1, 1)
    )
    result_no_valid <- metaanalysis_fx(no_valid_or_data)
    RUnit::checkTrue(nrow(result_no_valid) == 0)
}

# Run all the tests
testsuite <- defineTestSuite("TFEA.ChIP Tests", dirs = ".", testFileRegexp = "^test_.*\\.R$")
testResults <- runTestSuite(testsuite)

# Print the results
printTextProtocol(testResults)

