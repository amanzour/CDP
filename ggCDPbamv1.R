# CDP PLOT Version 1
# AUTHOR: AMIR MANZOUR

ggCDPbamv1 <- function(## INPUTS
  minTPM = 5,
  maxTPM = Inf,
  minintron = 5,
  minexon = 5,
  minExtensionRatio = 0.1,
  absminTPM = 1,
  absminintron = 1,
  absminexon = 1,
  absminExtensionRatio = 0.01,
  minTxSize = 500,
  minIntronSize = 200,
  minExonSize = 500,
  Extension = 1000,
  overlapcoef=0.95,
  RNAseqLengthNormalize=TRUE,
  PARCLIPLengthNormalize = TRUE,
  readlength=50,
  minoverlapOnIntrons=25,
  minoverlapOnExons=25,
  gtfpath = NULL,
  myWTpaths = NULL,
  myTreatmentpaths = NULL,
  myPARCLIPpath = NULL,
  returnAll = FALSE,
  MasterTableOriginal = NULL,
  MasterTable = NULL,
  ispairedendreadtext = "single-read",
  ispairedendread=FALSE,
  Only24chromosomes=TRUE,
  ReverseBamFileStrand=TRUE,
  FilterPseudogenes=TRUE,
  FilterAntisense=TRUE,
  TargetTolerance=50,
  ignoreParclipStrand = FALSE,
  #====CDF parameters
  includeBin0 = TRUE,
  plotByColumn="log2fC",
  binColumn="PARCLIPgene",
  histBreaks=1000,
  seehistogram=FALSE,
  binNames=c("bin0","bin1","bin2","bin3","bin4"),
  exprNormalization = TRUE,
  ## GGPLOT PARAMETERS
  fontFamily = "Helvetica",
  legendTitle = c("Bin","genes","median","pval"),
  legendTextSize = 5,
  legendPosition = c(1.5,0.5),
  limitsX = c(-2.5,2.5),
  pointSize = 1,
  lineSize = 0.2,
  panelBorderSize = 0.5,
  legendSpacing = 0.5,
  xlab = "Log2 Fold Change [Treatment/WT]",
  ylab = "Cumulative Distribution",
  axisTitleSize = 8,
  axisTextSize = 8,
  ggTitleSize = 5,
  AnalysisName = NULL
){
  
  ## //LIBRARIES
  ## ............
  
  suppressWarnings(suppressMessages(library(data.table)))
  suppressWarnings(suppressMessages(library(plyr)))
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(reshape2)))
  suppressWarnings(suppressMessages(library(GenomicFeatures)))
  suppressWarnings(suppressMessages(library(GenomicRanges)))
  suppressWarnings(suppressMessages(library(GenomeInfoDb)))
  suppressWarnings(suppressMessages(library(Rsamtools))) 
  suppressWarnings(suppressMessages(library(GenomicAlignments))) 
  suppressWarnings(suppressMessages(library(ggplot2)))
  suppressWarnings(suppressMessages(library(gplots)))
  suppressWarnings(suppressMessages(library(RColorBrewer)))
  suppressWarnings(suppressMessages(library(extrafont)))
  suppressWarnings(suppressMessages(library(grid)))
  suppressWarnings(suppressMessages(library(gridExtra)))
  suppressWarnings(suppressMessages(library(reporttools)))
  suppressWarnings(suppressMessages(library(rtracklayer)))
  # suppressWarnings(suppressMessages(library(MVN)))
  suppressWarnings(suppressMessages(library(ggrepel)))
  suppressWarnings(suppressMessages(library(plyr)))
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(ggpubr)))
  suppressWarnings(suppressMessages(library(plotly)))
  suppressWarnings(suppressMessages(library(idr)))
  suppressWarnings(suppressMessages(library(htmlwidgets)))
  suppressWarnings(suppressMessages(library(ggfortify)))
  suppressWarnings(suppressMessages(library(mclust)))
  # suppressWarnings(suppressMessages(library(DESeq2)))
  # Make Master Table
  if(is.null(MasterTable) == TRUE) {
    
    # make sure necessary files are in input
    if(is.null(myPARCLIPpath) == TRUE){stop("You forgot to include the PAR-CLIP XL bam file(s).")}
    if(is.null(myWTpaths) == TRUE){stop("You forgot to include WT RNA-seq bam files.")}
    if(is.null(myTreatmentpaths) == TRUE){stop("You forgot to include Treatment RNA-seq bam file(s).")}
    if(is.null(gtfpath) == TRUE){stop("You forgot to include the gtf file.")}
    
    message("constructing gene annotations...")
    
    # Reading GTF file
    gtftable <- rtracklayer::import(gtfpath)
    gtf=as.data.frame(gtftable)
    setDT(gtf)
    if ("gene_type" %in% colnames(gtf) & dim(gtf[type == "gene",c("gene_id","gene_name")])[1]>0){
    gene_names_types <- gtf[type == "gene",c("gene_id","gene_name","gene_type")]
    } else {
      gene_names_types <- gtf[type == "start_codon",c("gene_id","gene_name")]
      gene_names_types[["gene_type"]] <- gene_names_types[["gene_name"]]
    }
    
    setDT(gene_names_types)
    
    txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtfpath, format = "gtf", organism = "Homo sapiens")
    # PREPARE GENE COORDINATES
    genes_GRangelist <- GenomicFeatures::genes(txdb)
    
    if (Only24chromosomes){
      chromosomecount=24
    } else {
      chromosomecount=length(GenomeInfoDb::seqlevels(genes_GRangelist))
    }
    
    genes_GRangelist <- genes_GRangelist[GenomicRanges::seqnames(genes_GRangelist) %in% GenomeInfoDb::seqlevels(genes_GRangelist)[1:chromosomecount], ]
    
    transcriptInfo <- transcriptLengths(txdb = txdb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE)
    transcriptInfo <- as.data.frame(transcriptInfo)
    setDT(transcriptInfo)
    geneID <- transcriptInfo[, .SD[which.max(tx_len)], by = gene_id]
    geneID <- as.data.frame(geneID)
    setDT(geneID)
    
    geneInfo <-base::merge(geneID,gene_names_types, by = "gene_id", all.x = TRUE)
    geneInfo <-geneInfo[,c("gene_id","tx_len","gene_name","gene_type"),with=FALSE]
    m <- match(mcols(genes_GRangelist)$gene_id, geneInfo$gene_id)
    mcols(genes_GRangelist) <- cbind(mcols(genes_GRangelist), geneInfo[m, -1L, drop=FALSE])
    
    # INTRON COORDINATES
    myintrons <- GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE)
    # constrain introns on gene coordinates:
    myintrons_unlist <- unlist(myintrons)
    
    mygeneoverlaps <- findOverlaps(myintrons_unlist, genes_GRangelist)
    myintrons_unlist_ingenes <- myintrons_unlist[queryHits(mygeneoverlaps)]
    mcols(myintrons_unlist_ingenes)$gene_id <- genes_GRangelist[subjectHits(mygeneoverlaps)]$gene_id
    
    myintrons_unlist_ingenes_reduced <- reduce(myintrons_unlist_ingenes, drop.empty.ranges=TRUE,with.revmap=TRUE)
    myrevmap1 <- mcols(myintrons_unlist_ingenes_reduced)$revmap
    myrevmap1gene_id <- relist(mcols(myintrons_unlist_ingenes)[unlist(myrevmap1), ], myrevmap1)
    mcols(myintrons_unlist_ingenes_reduced)$gene_id <- unlist(lapply(myrevmap1gene_id, function(x) x[1]))
    
    # exclude exons
    exonsdub <- GenomicFeatures::exons(txdb, use.names=TRUE)
    exons_reduced <-reduce(exonsdub)
    myintrons_unlist_ingenes_reduced_diff <- setdiff(myintrons_unlist_ingenes_reduced, exons_reduced,ignore.strand=F)
    myintrons_unlist_ingenes_reduced_exoncleaned <- myintrons_unlist_ingenes_reduced[subjectHits(findOverlaps(myintrons_unlist_ingenes_reduced_diff, myintrons_unlist_ingenes_reduced))]
    
    myintrons_unlist_ingenes_reduced_exoncleaned_reduced <- reduce(myintrons_unlist_ingenes_reduced_exoncleaned, drop.empty.ranges=TRUE,with.revmap=TRUE)
    
    myrevmap2 <- mcols(myintrons_unlist_ingenes_reduced_exoncleaned_reduced)$revmap
    myrevmap2gene_id <- relist(mcols(myintrons_unlist_ingenes_reduced_exoncleaned)[unlist(myrevmap2), "gene_id"], myrevmap2)
    mcols(myintrons_unlist_ingenes_reduced_exoncleaned_reduced)$gene_id <- unlist(lapply(myrevmap2gene_id, function(x) x[1]))
    intronsfinalall <- myintrons_unlist_ingenes_reduced_exoncleaned_reduced
    
    #CDS coordinates
    cds <- GenomicFeatures::cdsBy(txdb, by="gene", use.names=FALSE)
    cds_unlist <- unlist(cds, use.names=TRUE)
    mcols(cds_unlist)$gene_id <- names(cds_unlist)
    
    #Exon coordinates
    exonsdub_unlist <- exonsdub
    mygeneoverlapsexon <- findOverlaps(exonsdub_unlist, genes_GRangelist)
    exonsdub_unlist_ingenes <- exonsdub_unlist[queryHits(mygeneoverlapsexon)]
    mcols(exonsdub_unlist_ingenes)$gene_id <- genes_GRangelist[subjectHits(mygeneoverlapsexon)]$gene_id
    
    exonsdub_unlist_ingenes_reduced <- reduce(exonsdub_unlist_ingenes, drop.empty.ranges=TRUE,with.revmap=TRUE)
    myrevmap1ex <- mcols(exonsdub_unlist_ingenes_reduced)$revmap
    myrevmap1exgene_id <- relist(mcols(exonsdub_unlist_ingenes)[unlist(myrevmap1ex), "gene_id"], myrevmap1ex)
    mcols(exonsdub_unlist_ingenes_reduced)$gene_id <- unlist(lapply(myrevmap1exgene_id, function(x) x[1]))
    exonsall <- exonsdub_unlist_ingenes_reduced
    
    #genes_GRangelist <- genes_GRangelist[GenomicRanges::seqnames(genes_GRangelist) %in% GenomeInfoDb::seqlevels(genes_GRangelist)[1:chromosomecount], ]
    # 3UTR extension
    fragment = readlength
    utr3_by_gene_flankraw1 <- flank(genes_GRangelist, width=fragment, start=FALSE, both=FALSE, use.names=TRUE, ignore.strand=FALSE)
    utr3_by_gene_flankraw <- flank(utr3_by_gene_flankraw1, width=Extension, start=FALSE, both=FALSE, use.names=TRUE, ignore.strand=FALSE)
    utr3_by_gene_flank <- utr3_by_gene_flankraw[-queryHits(GenomicRanges::findOverlaps(utr3_by_gene_flankraw, genes_GRangelist, type="any")),]
    utr3_by_gene_flank <- utr3_by_gene_flank[,c("gene_id")]
    
    utr3_by_tx <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names=FALSE)
    all_utr3 <- unlist(utr3_by_tx, use.names=FALSE)
    mcols(all_utr3)$tx_id <- rep(as.integer(names(utr3_by_tx)), lengths(utr3_by_tx))
    tx2gene <- mcols(transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id")))
    tx2gene$gene_id <- as.character(tx2gene$gene_id)
    m <- match(mcols(all_utr3)$tx_id, tx2gene$tx_id)
    mcols(all_utr3) <- cbind(mcols(all_utr3), tx2gene[m, -1L, drop=FALSE])
    m <- match(unique(all_utr3),all_utr3)
    all_utr3_unique <- all_utr3[m]
    
    utr5_by_tx <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names=FALSE)
    all_utr5 <- unlist(utr5_by_tx, use.names=FALSE)
    mcols(all_utr5)$tx_id <- rep(as.integer(names(utr5_by_tx)), lengths(utr5_by_tx))
    tx2gene <- mcols(transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id")))
    tx2gene$gene_id <- as.character(tx2gene$gene_id)
    m <- match(mcols(all_utr5)$tx_id, tx2gene$tx_id)
    mcols(all_utr5) <- cbind(mcols(all_utr5), tx2gene[m, -1L, drop=FALSE])
    m <- match(unique(all_utr5),all_utr5)
    all_utr5_unique <- all_utr5[m]
    
    
    # Reading all bam files
    message("reading WT RNA-seq bam files... ")
    
    if (ispairedendread){
      ispairedendreadtext = "paired-end-read"
    }
    myWTbams <- lapply(myWTpaths, function (i){
      if (ispairedendread){
        myparam <- ScanBamParam(flag=scanBamFlag(isProperPair=TRUE,
                                                 isDuplicate=NA,
                                                 isSecondaryAlignment=FALSE))
        
        retvalpairs <- GenomicAlignments::readGAlignmentPairs(i, 
                                                              use.names = TRUE, 
                                                              strandMode = 1, 
                                                              param = myparam)
        retval <- GenomicAlignments::granges(retvalpairs) 
        
        if (ReverseBamFileStrand){
          retval <- invertStrand(retval[GenomicRanges::seqnames(retval) %in% GenomeInfoDb::seqlevels(genes_GRangelist)[1:chromosomecount], ])
        } else {
          retval <- retval[GenomicRanges::seqnames(retval) %in% GenomeInfoDb::seqlevels(genes_GRangelist)[1:chromosomecount], ]
        }
      } else {
        myparam <- ScanBamParam(flag=scanBamFlag(isDuplicate=NA,
                                                 isSecondaryAlignment=FALSE))
        
        if (ReverseBamFileStrand){
          retval <- invertStrand(GenomicAlignments::readGAlignments(i, 
                                                                    use.names = FALSE, 
                                                                    param = myparam))
        } else {
          retval <- GenomicAlignments::readGAlignments(i, 
                                                       use.names = FALSE, 
                                                       param = myparam)
        }
      }
      message(paste0("...imported ",ispairedendreadtext," WT bam file ",i))
      return(retval)
    })
    message("OK!")
    message("reading Treatment RNA-seq bam files... ")
    myTreatmentbams <- lapply(myTreatmentpaths, function (i){
      if (ispairedendread){
        myparam <- ScanBamParam(flag=scanBamFlag(isProperPair=TRUE,
                                                 isDuplicate=NA,
                                                 isSecondaryAlignment=FALSE))
        
        retvalpairs <- GenomicAlignments::readGAlignmentPairs(i, 
                                                              use.names = TRUE, 
                                                              strandMode = 1, 
                                                              param = myparam)
        retval <- GenomicAlignments::granges(retvalpairs) 
        retval <- invertStrand(retval[GenomicRanges::seqnames(retval) %in% GenomeInfoDb::seqlevels(genes_GRangelist)[1:chromosomecount], ])
      } else {
        myparam <- ScanBamParam(flag=scanBamFlag(isDuplicate=NA,
                                                 isSecondaryAlignment=FALSE))
        
        retval <- invertStrand(GenomicAlignments::readGAlignments(i, 
                                                                  use.names = FALSE, 
                                                                  param = myparam))
      }
      message(paste0("...imported ",ispairedendreadtext," Treatment bam file ",i))
      return(retval)
    })
    
    message("OK!")
    message("reading the Parclip XL-reads bam file... ")
    # READE PARCLIP
    myPARCLIPbam <- GenomicAlignments::readGAlignments(myPARCLIPpath, use.names = TRUE)
    myPARCLIPbam <- myPARCLIPbam[GenomicRanges::seqnames(myPARCLIPbam) %in% GenomeInfoDb::seqlevels(genes_GRangelist)[1:chromosomecount], ]
    message("OK!")
    
    message("mapping RNA-seq reads to gene annotation regions...")
    #Gene
    #minoverlapOnGenes=round(overlapcoef*readlength)
    for (repID in 1:length(myWTbams)){
      mcols(genes_GRangelist)[[paste0("WTgene",repID)]] <- GenomicRanges::countOverlaps(genes_GRangelist, myWTbams[[repID]], type = "any", maxgap = TargetTolerance, ignore.strand = FALSE)
    }
    for (repID in 1:length(myTreatmentbams)){
      mcols(genes_GRangelist)[[paste0("Treatmentgene",repID)]] <- GenomicRanges::countOverlaps(genes_GRangelist, myTreatmentbams[[repID]], type = "any", maxgap = TargetTolerance, ignore.strand = FALSE)
    }
    mcols(genes_GRangelist)[["PARCLIPgene"]] <- GenomicRanges::countOverlaps(genes_GRangelist, myPARCLIPbam, type = "any", ignore.strand = ignoreParclipStrand)
    genes_df <- as.data.frame(genes_GRangelist)
    setDT(genes_df)
    
    #Exon
    for (repID in 1:length(myWTbams)){
      mcols(exonsall)[[paste0("WTexon",repID)]] <- GenomicRanges::countOverlaps(exonsall, myWTbams[[repID]], type = "any", minoverlap = minoverlapOnExons, ignore.strand = FALSE)
    }
    for (repID in 1:length(myTreatmentbams)){
      mcols(exonsall)[[paste0("Treatmentexon",repID)]] <- GenomicRanges::countOverlaps(exonsall, myTreatmentbams[[repID]], type = "any", minoverlap = minoverlapOnExons, ignore.strand = FALSE)
    }
    mcols(exonsall)[["PARCLIPexon"]] <- GenomicRanges::countOverlaps(exonsall, myPARCLIPbam, type = "any", ignore.strand = ignoreParclipStrand)
    
    exonsalldf <- as.data.frame(exonsall)
    setDT(exonsalldf)
    exons_gene <- exonsalldf[, lapply(.SD,sum), .SDcols = which(grepl("WTexon|Treatmentexon|PARCLIPexon|width",colnames(exonsalldf))), by = "gene_id"]
    colnames(exons_gene)[which(colnames(exons_gene)=="width")] <- "exonwidth"
    #Intron
    for (repID in 1:length(myWTbams)){
      mcols(intronsfinalall)[[paste0("WTintron",repID)]] <- GenomicRanges::countOverlaps(intronsfinalall, myWTbams[[repID]], type = "any", minoverlap = minoverlapOnIntrons, ignore.strand = FALSE)
    }
    for (repID in 1:length(myTreatmentbams)){
      mcols(intronsfinalall)[[paste0("Treatmentintron",repID)]] <- GenomicRanges::countOverlaps(intronsfinalall, myTreatmentbams[[repID]], type = "any", minoverlap = minoverlapOnIntrons, ignore.strand = FALSE)
    }
    mcols(intronsfinalall)[["PARCLIPintron"]] <- GenomicRanges::countOverlaps(intronsfinalall, myPARCLIPbam, type = "any", ignore.strand = ignoreParclipStrand)
    
    intronsfinalalldf <- as.data.frame(intronsfinalall)
    setDT(intronsfinalalldf)
    introns_gene <- intronsfinalalldf[, lapply(.SD,sum), .SDcols = which(grepl("WTintron|Treatmentintron|PARCLIPintron|width",colnames(intronsfinalalldf))), by = "gene_id"]
    colnames(introns_gene)[which(colnames(introns_gene)=="width")] <- "intronwidth"
    # Extension
    minoverlapOnExtensions=round(overlapcoef*readlength)
    for (repID in 1:length(myWTbams)){
      mcols(utr3_by_gene_flank)[[paste0("WTExtension",repID)]] <- GenomicRanges::countOverlaps(utr3_by_gene_flank, myWTbams[[repID]], type = "any", minoverlap = minoverlapOnExtensions, ignore.strand = FALSE)
    }
    for (repID in 1:length(myTreatmentbams)){
      mcols(utr3_by_gene_flank)[[paste0("TreatmentExtension",repID)]] <- GenomicRanges::countOverlaps(utr3_by_gene_flank, myTreatmentbams[[repID]], type = "any", minoverlap = minoverlapOnExtensions, ignore.strand = FALSE)
    }
    utr3_by_gene_flank_df <- as.data.frame(utr3_by_gene_flank)
    colnames(utr3_by_gene_flank_df)[which(colnames(utr3_by_gene_flank_df)=="width")] <- "Extensionwidth"
    setDT(utr3_by_gene_flank_df)
    
    #CDS
    mcols(cds_unlist)[["PARCLIPCDS"]] <- GenomicRanges::countOverlaps(cds_unlist, myPARCLIPbam, type = "any", ignore.strand = ignoreParclipStrand)
    cds_unlistdf <- as.data.frame(cds_unlist, row.names = NULL)
    setDT(cds_unlistdf)
    cds_gene <- cds_unlistdf[, lapply(.SD,sum), .SDcols = which(grepl("PARCLIPCDS",colnames(cds_unlistdf))), by = "gene_id"]
    
    #3UTR
    mcols(all_utr3_unique)[["PARCLIP3UTR"]] <- GenomicRanges::countOverlaps(all_utr3_unique, myPARCLIPbam, type = "any", maxgap = TargetTolerance, ignore.strand = ignoreParclipStrand)
    
    all_utr3_uniquedf <- as.data.frame(all_utr3_unique)
    setDT(all_utr3_uniquedf)
    utr3_gene <- all_utr3_uniquedf[, lapply(.SD,sum), .SDcols = which(grepl("PARCLIP3UTR",colnames(all_utr3_uniquedf))), by = "gene_id"]
    
    #5UTR
    mcols(all_utr5_unique)[["PARCLIP5UTR"]] <- GenomicRanges::countOverlaps(all_utr5_unique, myPARCLIPbam, type = "any", maxgap = TargetTolerance, ignore.strand = ignoreParclipStrand)
    
    all_utr5_uniquedf <- as.data.frame(all_utr5_unique)
    setDT(all_utr5_uniquedf)
    utr5_gene <- all_utr5_uniquedf[, lapply(.SD,sum), .SDcols = which(grepl("PARCLIP5UTR",colnames(all_utr5_uniquedf))), by = "gene_id"]
    
    # UPTO HERE NOW PARCLIP
    
    message("integrating information to a preliminary master table: MasterTableOriginal..")
    # Make MasterTable 
    MasterTable0 <- base::merge(genes_df,exons_gene, by = "gene_id", all.x = TRUE)
    MasterTable1 <- base::merge(MasterTable0,introns_gene, by = "gene_id", all.x = TRUE)
    MasterTable1cds <- base::merge(MasterTable1,cds_gene, by = "gene_id", all.x = TRUE)
    MasterTable2 <- base::merge(MasterTable1cds,utr3_gene, by = "gene_id", all.x = TRUE)
    MasterTable3 <- base::merge(MasterTable2,utr5_gene, by = "gene_id", all.x = TRUE)
    MasterTableOriginal <- base::merge(MasterTable3,utr3_by_gene_flank_df[,c("gene_id","Extensionwidth",colnames(utr3_by_gene_flank_df)[which(grepl("WTExtension|TreatmentExtension",colnames(utr3_by_gene_flank_df)))]),with=FALSE], by = "gene_id", all.x = TRUE)
    setDT(MasterTableOriginal)
    MasterTableOriginal
    #numericcolnames <- colnames(MasterTableOriginal)[-which(colnames(MasterTableOriginal) %in% c("gene_id","seqnames","start","end","width","strand","tx_len","gene_name","gene_type"))]
    #for(j in (numericcolnames)){data.table::set(MasterTableOriginal, i = which(is.na(MasterTableOriginal[[j]])), j = j, value = 0)}
    message("removing all pseudogenes and RNA genes from analysis...")
    
    MasterTableOriginal <- MasterTableOriginal[!grepl("^RPL|^RPS|^EEEF1A|^EEF1A|^LENG7",gene_name) & !grepl("RNA$",gene_type) ,]
    if (FilterPseudogenes) {
      MasterTableOriginal <- MasterTableOriginal[!grepl("pseudogene$",gene_type), ]
    }
    if (FilterAntisense){
      MasterTableOriginal <- MasterTableOriginal[!grepl("antisense$",gene_type), ]
    }
    message("OK!")
    
    MasterTable <- MasterTableOriginal
    
    message("averaging across ",length(myWTpaths)," WT replicates and ",length(myTreatmentpaths)," Treatment replicates...  computing gene TPM/RPKM, downstream extension coverage, exonic, and intronic reads... taking log2 differentials between WT and Treatment expressions... ")
    
    # GET EXTENSION TPM READY
    extensioncolnames <- colnames(MasterTable)[which(grepl("WTExtension|TreatmentExtension",colnames(MasterTable)))]
    for (k in extensioncolnames) {data.table::set(x=MasterTable, j=k, value = (10^6/sum(MasterTable[[gsub("Extension","gene",k)]],na.rm = TRUE))*MasterTable[[k]]/1)}
    
    # GET GENE TPM READY
    if (!RNAseqLengthNormalize){
      genecolnames <- colnames(MasterTable)[which(grepl("WTgene|Treatmentgene",colnames(MasterTable)))]
      for (k in genecolnames) {data.table::set(x=MasterTable, j=k, value = (10^6/sum(MasterTable[[k]],na.rm = TRUE))*MasterTable[[k]]/1)}
      #for(j in (genecolnames)){data.table::set(MasterTable, i = which(is.nan(MasterTable[[j]])), j = j, value = 0)}
      
      # GET EXON PM READY
      exoncolnames <- colnames(MasterTable)[which(grepl("WTexon|Treatmentexon",colnames(MasterTable)))]
      for (k in exoncolnames) {data.table::set(x=MasterTable, j=k, value = (10^6/sum(MasterTable[[k]],na.rm = TRUE))*MasterTable[[k]]/1)}
      
      # GET INTRON PM READY
      introncolnames <- colnames(MasterTable)[which(grepl("WTintron|Treatmentintron",colnames(MasterTable)))]
      for (k in introncolnames) {data.table::set(x=MasterTable, j=k, value = (10^6/sum(MasterTable[[k]],na.rm = TRUE))*MasterTable[[k]]/1)}
    } else {
      message("normalizing RNA-seq expression by gene/exon/intron lengths")
      genecolnames <- colnames(MasterTable)[which(grepl("WTgene|Treatmentgene",colnames(MasterTable)))]
      for (k in genecolnames) {data.table::set(x=MasterTable, j=k, value = (10^9/sum(MasterTable[[k]]/MasterTable[["width"]],na.rm = TRUE))*MasterTable[[k]]/MasterTable[["width"]])}
      #for(j in (genecolnames)){data.table::set(MasterTable, i = which(is.nan(MasterTable[[j]])), j = j, value = 0)}
      
      # GET EXON PM READY
      exoncolnames <- colnames(MasterTable)[which(grepl("WTexon|Treatmentexon",colnames(MasterTable)))]
      for (k in exoncolnames) {data.table::set(x=MasterTable, j=k, value = (10^9/sum(MasterTable[[k]]/MasterTable[["exonwidth"]],na.rm = TRUE))*MasterTable[[k]]/MasterTable[["exonwidth"]])}
      
      # GET INTRON PM READY
      introncolnames <- colnames(MasterTable)[which(grepl("WTintron|Treatmentintron",colnames(MasterTable)))]
      for (k in introncolnames) {data.table::set(x=MasterTable, j=k, value = (10^9/sum(MasterTable[[k]]/MasterTable[["intronwidth"]],na.rm = TRUE))*MasterTable[[k]]/MasterTable[["intronwidth"]])}
    }
    
    # AVERAGING TPM
    WTgenecolnames <- colnames(MasterTable)[which(grepl("WTgene",colnames(MasterTable)))]
    MasterTable[, WTTPM := rowMeans(MasterTable[,WTgenecolnames,with=FALSE],na.rm = TRUE)]
    Treatmentgenecolnames <- colnames(MasterTable)[which(grepl("Treatmentgene",colnames(MasterTable)))]
    MasterTable[, TreatmentTPM := rowMeans(MasterTable[,Treatmentgenecolnames,with=FALSE],na.rm = TRUE)]
    
    #AVERAGING EXTENSION TPM
    WTextensioncolnames <- colnames(MasterTable)[which(grepl("WTExtension",colnames(MasterTable)))]
    MasterTable[, WTExtensionTPM := rowMeans(MasterTable[,WTextensioncolnames,with=FALSE],na.rm = TRUE)]
    Treatmentextensioncolnames <- colnames(MasterTable)[which(grepl("TreatmentExtension",colnames(MasterTable)))]
    MasterTable[, TreatmentExtensionTPM := rowMeans(MasterTable[,Treatmentextensioncolnames,with=FALSE],na.rm = TRUE)]
    # RATIO OF EXTENSION TO GENE TPM
    MasterTable[, WTExtensionRatio := (MasterTable[["WTExtensionTPM"]]/MasterTable[["Extensionwidth"]])/(MasterTable[["WTTPM"]]/MasterTable[["width"]])]
    MasterTable[, TreatmentExtensionRatio := (MasterTable[["TreatmentExtensionTPM"]]/MasterTable[["Extensionwidth"]])/(MasterTable[["TreatmentTPM"]]/MasterTable[["width"]])]
    #=========================
    #AVERAGING EXON
    WTexoncolnames <- colnames(MasterTable)[which(grepl("WTexon",colnames(MasterTable)))]
    MasterTable[, WTexonmean := rowMeans(MasterTable[,WTexoncolnames,with=FALSE],na.rm = TRUE)]
    Treatmentexoncolnames <- colnames(MasterTable)[which(grepl("Treatmentexon",colnames(MasterTable)))]
    MasterTable[, Treatmentexonmean := rowMeans(MasterTable[,Treatmentexoncolnames,with=FALSE],na.rm = TRUE)]
    #AVERAGING INTRON
    WTintroncolnames <- colnames(MasterTable)[which(grepl("WTintron",colnames(MasterTable)))]
    MasterTable[, WTintronmean := rowMeans(MasterTable[,WTintroncolnames,with=FALSE],na.rm = TRUE)]
    Treatmentintroncolnames <- colnames(MasterTable)[which(grepl("Treatmentintron",colnames(MasterTable)))]
    MasterTable[, Treatmentintronmean := rowMeans(MasterTable[,Treatmentintroncolnames,with=FALSE],na.rm = TRUE)]
    #=========================
    MasterTableReplicates <- MasterTable
    MasterTable <- MasterTable[,c("seqnames","start","end","gene_id","gene_name","strand","gene_type","tx_len","width","exonwidth","intronwidth","Extensionwidth","WTTPM","TreatmentTPM","WTexonmean","Treatmentexonmean","WTintronmean","Treatmentintronmean","WTExtensionTPM","TreatmentExtensionTPM","WTExtensionRatio","TreatmentExtensionRatio","PARCLIPgene","PARCLIP5UTR","PARCLIPCDS","PARCLIPintron","PARCLIP3UTR")]
    
    MasterTable[, log2fC := log2((MasterTable[["TreatmentTPM"]]+absminTPM)/(MasterTable[["WTTPM"]]+absminTPM))]
    MasterTable[, log2fCexon := log2((MasterTable[["Treatmentexonmean"]]+absminexon)/(MasterTable[["WTexonmean"]]+absminexon))]
    MasterTable[, log2fCintron := log2((MasterTable[["Treatmentintronmean"]]+absminintron)/(MasterTable[["WTintronmean"]]+absminintron))]
    MasterTable[, log2fCextension := log2((MasterTable[["TreatmentExtensionRatio"]]+absminExtensionRatio)/(MasterTable[["WTExtensionRatio"]]+absminExtensionRatio))]
    
    message("DONE! ... with part 1")
    
    
    # return Master Table
    return(list(MasterTable, MasterTableReplicates))
  } else {
    
    ## //INPUT CHECKING
    ## .................
    
    
    # filter, subset and extract the columns for plotting and calculations
    MasterTable <- MasterTable[WTTPM >= absminTPM & 
                                 TreatmentTPM >= absminTPM & 
                                 tx_len >= minTxSize &
                                 (WTTPM <= maxTPM & TreatmentTPM <= maxTPM) &
                                 (WTTPM >= minTPM | TreatmentTPM >= minTPM) ]
    message(paste0("filtered MasterTable for transcript length of at least ",minTxSize,"nt."))
    message(paste0("filtered MasterTable for RPKM of at least ",absminTPM," in both samples and at least ",minTPM," in either sample."))
    message(paste0("filtered MasterTable for RPKM of at most ",maxTPM," in both samples."))
    
    if (plotByColumn == "log2fCexon"){
      MasterTable <- MasterTable[WTexonmean >= absminexon & 
                                   Treatmentexonmean >= absminexon & 
                                   exonwidth >= minExonSize &
                                   (WTexonmean >= minexon | Treatmentexonmean >= minexon) ]
      message(paste0("additionally filtered MasterTable for exon length of at least ",minExonSize,"nt."))
      message(paste0("additionally filtered MasterTable for average exonic RPKM of at least ",absminexon, " in both samples and at least ",minintron," in either sample."))
    }
    
    if (plotByColumn == "log2fCintron"){
      MasterTable <- MasterTable[WTintronmean >= absminintron & 
                                   Treatmentintronmean >= absminintron & 
                                   intronwidth >= minIntronSize &
                                   (WTintronmean >= minintron | Treatmentintronmean >= minintron) ]
      message(paste0("additionally filtered MasterTable for intron length of at least ",minIntronSize,"nt."))
      message(paste0("additionally filtered MasterTable for average intronic RPKM of at least ",absminintron, " in both samples and at least ",minintron," in either sample."))
    }
    
    if (plotByColumn == "log2fCextension"){
      MasterTable <- MasterTable[WTExtensionRatio >= absminExtensionRatio & TreatmentExtensionRatio >= absminExtensionRatio &
                                   (WTExtensionRatio >= minExtensionRatio | TreatmentExtensionRatio >= minExtensionRatio)]
      message(paste0("additionally filtered MasterTable for downstream expression of at least ",absminExtensionRatio," of corresponding gene expression in both samples and at least ",minExtensionRatio," in either sample."))
    }
    # (WTExtensionRatio > absminExtensionRatio & TreatmentExtensionRatio > absminExtensionRatio) &
    # & (WTExtensionRatio > minExtensionRatio | TreatmentExtensionRatio > minExtensionRatio)
    # write.csv(MasterTable,"MasterTable.csv",sep = "\t")
    
    # RNA SEQ LENGTH NORMALIZATION
    
    # PARCLIP BINNING
    #in the next version, normalize the binning by WTTPM, WTexonmean, and WTintronmean
    if (exprNormalization){
      if (!PARCLIPLengthNormalize){
        message(paste0("gene binning of: ",plotByColumn," for target region: ",binColumn,", using log2 [ (10^6*XL/sum(XL) / (Gene Expression in WT) ] gene distribution..."))    
        MasterTable[, target := log2(10^6*(MasterTable[[binColumn]]/sum(MasterTable[[binColumn]],na.rm = TRUE))/(MasterTable[["WTTPM"]]+0))]
      } else {
        message(paste0("gene binning of: ",plotByColumn," for target region: ",binColumn,", using log2 [ (10^9*XL/sum(XL) / (Length-Normalized Gene Expression  in WT) ] gene distribution..."))    
        MasterTable[, target := log2(10^9*((MasterTable[[binColumn]]/MasterTable[["width"]])/sum(MasterTable[[binColumn]]/MasterTable[["width"]],na.rm = TRUE))/(MasterTable[["WTTPM"]]+0))]
      }
    } else {
      if (!PARCLIPLengthNormalize){
        message(paste0("gene binning of: ",plotByColumn," for target region: ",binColumn,", using log2 [ (10^6*XL/sum(XL) ] gene distribution..."))    
        MasterTable[, target := log2(10^6*(MasterTable[[binColumn]]/sum(MasterTable[[binColumn]],na.rm = TRUE)))]
      } else {
        message(paste0("gene binning of: ",plotByColumn," for target region: ",binColumn,", using log2 [ (10^9*(XL/GeneLength)/sum((XL/GeneLength)) ] gene distribution..."))    
        MasterTable[, target := log2(10^9*((MasterTable[[binColumn]]/MasterTable[["width"]])/sum(MasterTable[[binColumn]]/MasterTable[["width"]],na.rm = TRUE))/(1+0))]
      }
    }
    
    
    targetmean <- mean(unlist(MasterTable[is.finite(MasterTable[["target"]]),"target",with=FALSE]),na.rm = TRUE)
    targetsd <- sd(unlist(MasterTable[is.finite(MasterTable[["target"]]),"target",with=FALSE]),na.rm = TRUE)
    q1 <- targetmean + targetsd
    q2 <- targetmean
    q3 <- targetmean - targetsd
    MasterTable[, binning := binNames[1]]
    MasterTable[is.finite(target),"binning"] <- binNames[2]
    MasterTable[target > q3 & target <= q2 ,"binning"] <- binNames[3]
    MasterTable[target > q2 & target <= q1 ,"binning"] <- binNames[4]
    MasterTable[target > q1,"binning"] <- binNames[5]
    
    message("OK!")
    
    if (!includeBin0){
      MasterTable <- MasterTable[binning != "bin0",]
    }
    # subset the master table
    MT <- MasterTable[, c("gene_name","gene_type","log2fC","log2fCexon","log2fCintron","log2fCextension","target","binning"), with=FALSE]
    
    
    # helper histogram
    if (seehistogram == TRUE) {
      print(hist(unlist(MT[is.finite(MT[["target"]]),"target",with=FALSE]), 
                 breaks = histBreaks,main=paste0(binColumn," target distribution on genes")))
    }
    
    
    ## //CREATE SUMMARY DATA
    ## ......................
    
    plot_colors <- c("black", rev(brewer.pal(length(binNames[-1]), "Spectral")))
    ggtitle = paste(AnalysisName,"\nCDP plotted by: ",plotByColumn," binned by: ", binColumn, sep = "")
    
    # select the data frame
    DF <- data.frame(MT[,c(plotByColumn,"target","binning"), with=FALSE])
    
    
    # calculate medians
    DFmed <- data.frame(bin = binNames,
                        x = unlist(lapply(binNames,function(i){signif(median(DF[DF[,"binning"] == i,][,plotByColumn]),digits = 2)})),
                        y = rep(0.01,length(binNames)))
    
    message("...generated median table for all the bins.")
    
    # calculate legend
    binRange <- lapply(binNames, function(i){
      imin <- signif(min(DF[DF[,"binning"] == i,][,"target"],na.rm = TRUE),2)
      imax <- signif(max(DF[DF[,"binning"] == i,][,"target"],na.rm = TRUE),2)
      return(paste(imin," , ",imax,sep = ""))})
    
    DFks <-  data.frame(bin = binNames, 
                        binRange = unlist(binRange),
                        n = as.character(unlist(lapply(binNames, function(i){nrow(DF[DF[,"binning"] == i,])}))),
                        m = DFmed$x,
                        x = c(unlist(lapply(binNames,function(i){retval=NA;if ((!includeBin0) | (dim(DF[DF[,"binning"] == binNames[1],])[1]==0)) {return(NA)} else {retval=format.pval(
                          ks.test(DF[DF[,"binning"] == binNames[1],][,plotByColumn],DF[DF[,"binning"] == i,][,plotByColumn],
                                  alternative = "two.sided")$p.val,digits = 3);return(retval)}}))))
    
    
    message("...calculated p-values using ks.test")
    
    ## //PLOTTING THE RESULTS
    ## .......................
    
    
    p<-ggplot(DF, aes_string(x = plotByColumn)) + 
      stat_ecdf(aes(colour =binning), geom = "step", size = lineSize) + 
      xlab(label = xlab) + 
      ylab(label = ylab) + 
      ggtitle(label = ggtitle) +
      theme_bw() + 
      scale_x_continuous(limits = limitsX) +
      scale_y_continuous(expand = c(0.01,0)) +
      
      scale_color_manual(values = plot_colors, name = paste(legendTitle[1],"; ",legendTitle[2],"; ",legendTitle[3],"; ",legendTitle[4],sep = ""),
                         labels = paste(DFks[,2],"; ",DFks[,3],"; ",DFks[,4],"; ", DFks[,5], sep = "")) +
      
      geom_point(data = DFmed, aes(x = x, y = y, colour = factor(bin)), size = pointSize) +
      
      theme(plot.title = element_text(size = ggTitleSize, family = fontFamily, hjust = 0.5),
            panel.grid = element_blank(), 
            panel.border = element_rect(colour = "black", size = panelBorderSize),
            legend.position = legendPosition,
            legend.title = element_text(size = legendTextSize),
            legend.background = element_blank(),
            legend.key.size = unit(legendSpacing, 'lines'),
            legend.text = element_text(size = legendTextSize, family = fontFamily),
            axis.title.x = element_text(size = axisTitleSize, family = fontFamily,vjust = -1),
            axis.text.x = element_text(size = axisTextSize, family = fontFamily, face = "plain"),
            axis.title.y = element_text(size = axisTitleSize, family = fontFamily,vjust = 2),
            axis.text.y = element_text(size = axisTextSize, family = fontFamily, face = "plain"),
            plot.margin = margin(2,2,2,2,"cm"),
            aspect.ratio = 1)
    # plot.margin = margin(1,0.25,1,0.25,"cm")
    if (returnAll == TRUE) {
      dataList <- list(MasterTable, MT, DFmed, DFks, p, plot_colors)
      names(dataList) <- c("AveragedValues", "binnedTable", "MedianTable", "Kolmogorov-SmirnovTable","ggPlot","plotColors")
      return(dataList)
    } else {
      return(p)
    }
  }
  message("DONE!")
}