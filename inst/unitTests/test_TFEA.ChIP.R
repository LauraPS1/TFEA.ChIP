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
    RUnit::checkTrue(
        is.na(GeneID2entrez("potato",return.Matrix = TRUE)$ENTREZ.ID))
}

test_preprocessInputData<-function(){
    data("hypoxia_DESeq","hypoxia",package = "TFEA.ChIP")
    # Checking output size
    RUnit::checkEquals(
        nrow(preprocessInputData(hypoxia)),
        length(GeneID2entrez(hypoxia$Genes)))

    RUnit::checkEquals(
        nrow(preprocessInputData(hypoxia_DESeq)),
        length(GeneID2entrez(rownames(hypoxia_DESeq))))

    # Cheking every gene in output has been properly translated
    RUnit::checkTrue(!(FALSE %in% (
        preprocessInputData(hypoxia)$Genes %in% GeneID2entrez(hypoxia$Genes))))
    RUnit::checkTrue(!(FALSE %in% (
        preprocessInputData(hypoxia_DESeq)$Genes %in%
        GeneID2entrez(rownames(hypoxia_DESeq)))))

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
        nrow(MetaData[grepl("wg",MetaData$Accession),]))
    # TF "JUND" only
    RUnit::checkEquals(
        nrow(get_chip_index(TFfilter = "JUND")),
        nrow(MetaData[MetaData$TF=="JUND",]))
}

test_contingency_matrix<-function(){
    data("Genes.Upreg",package = "TFEA.ChIP")

    # Sum of every contingency matrix. Since no control gene list
    # is provided, all genes not included in the test list will
    # be used as control.
    cont_mat<-contingency_matrix(Genes.Upreg)
    RUnit::checkTrue(!(FALSE %in%
        (sapply(cont_mat,sum) == 23056))
    )
    # Size of the output when selecting ChIP-Seq datasets
    chip_index<-get_chip_index(TFfilter = "EPAS1")
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
    data("CM_list",package = "TFEA.ChIP")

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
        length(ARNT.gr<-txt2GR(ARNT.peaks.bed,"macs",ARNT.metadata)),
        nrow(ARNT.peaks.bed[ARNT.peaks.bed$X.10.log10.pvalue.<50]))
    # Checking that missing input variables produces errors
    RUnit::checkException(
        ARNT.gr<-txt2GR(ARNT.peaks.bed,fileMetaData = ARNT.metadata))
    RUnit::checkException(
        ARNT.gr<-txt2GR(ARNT.peaks.bed,"macs"))
    # Checking incorrect input variables produces errors
    RUnit::checkException(
        ARNT.gr<-txt2GR(ARNT.peaks.bed,"macs",ARNT.metadata[,1:3]))

}

test_GR2tfbs_db<-function(){
    data("DnaseHS_db","gr.list", package="TFEA.ChIP")
    # Using this toy datasets, the expected result is a
    # list of one gene set with only two gene IDs: 2782 and 23261
    RUnit::checkEquals(length(GR2tfbs_db(DnaseHS_db, gr.list)),1)
    RUnit::checkEquals(
        GR2tfbs_db(DnaseHS_db, gr.list)[[1]]@geneIds,
        c("2782","23261"))
}

test_GSEA_run<-function(){
    data("Entrez.gene.IDs",package = "TFEA.ChIP")
    # Selecting a small number of ChIP-Seq datasets to save time
    chip_index<-get_chip_index(TFfilter = c("EPAS1","ARNT"))
    # Checking enrichment table
    RUnit::checkEquals(
        GSEA_run(Entrez.gene.IDs,chip_index)$Accession,
        chip_index$Accession)
    expected_ES<-c(0.479650,0.449680,0.232260,0.298820,0.346470,
        0.226810,0.226920,0.404960,0.172440,0.203680, -0.047051,
        0.159360,0.266290,0.216330,-0.050501)
    RUnit::checkEquals(
        GSEA_run(Entrez.gene.IDs,chip_index)$ES,
        expected_ES)
    # Checking output with RES and indicator
    RUnit::checkTrue(is.list(
        GSEA_run(Entrez.gene.IDs,chip_index[1:2,],get.RES = TRUE)))


}

test_makeTFBSmatrix<-function(){
    data("tfbs.database","Entrez.gene.IDs",package = "TFEA.ChIP")
    # Checking output size
    RUnit::checkEquals(
        nrow(makeTFBSmatrix(Entrez.gene.IDs,tfbs.database)),
        length(Entrez.gene.IDs))
    RUnit::checkEquals(
        ncol(makeTFBSmatrix(Entrez.gene.IDs,tfbs.database)),
        length(tfbs.database))
    # Checking at least one gene is assigned to every ChIP-Seq
    RUnit::checkTrue(!(FALSE %in% (
        colSums(makeTFBSmatrix(Entrez.gene.IDs,tfbs.database))>0)))
    # Checking gene assignment is correct
    tmp<-makeTFBSmatrix(Entrez.gene.IDs,tfbs.database)[,1]
    RUnit::checkTrue(!(FALSE %in% (
        names(tmp[tmp==1]) %in% tfbs.database[[1]]@geneIds)))
}


















