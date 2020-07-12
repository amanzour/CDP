# CDP PLOT Version 1
# AUTHOR: AMIR MANZOUR
# https://bookdown.org/yihui/rmarkdown/parameterized-reports.html
rm(list=ls(all=TRUE))
#.rs.restartR() # BEST WAY TO CLEAR MEMORY
gc()
RunOnServer = FALSE  # TRUE: when running NIAMS server. FALSE: when running on machine


# setwd(rootpath)

# path to the program. you can keep this as it is and only change paramters below.
#source(paste0(rootpath,"bin/ggCDPbamv1.R"))
source(paste0("ggCDPbamv1.R"))
#debugSource(paste0("ggCDPbamv1.R"))
#debugSource(paste0(rootpath,"bin/ggCDPbamv1.R"))
myaccount="/home/manzouro/"
if (RunOnServer) {
  rootpath=myaccount
} else {
  rootpath="/Volumes/hpchome/"
  rootpath="/Users/manzouro/Documents/Research/StandardSoftware/"
}
myPARCLIPpath = paste0(rootpath,"HNRNPK_P_Cyt_224_226_RPI36.aligned_TtoC.sorted.bam")
# myPARCLIPpath=paste0(rootpath,gsub(myaccount,"",myPARCLIPpath))

myWTpaths = c("AF_P_mCh_V_CKDL200153406-1a-DY0088-AK1682_HC3L5BBXX_L7Aligned.sortedByCoord.out.bam")
myWTpaths = paste0(rootpath,myWTpaths)
myTreatmentpaths = c("AF_KO1_mCh_V_CKDL200153406-1a-DY0088-AK1544_HC3L5BBXX_L7Aligned.sortedByCoord.out.bam")
myTreatmentpaths = paste0(rootpath,myTreatmentpaths)
mygtfpath=c("genes_no_mir_snor_hist.gtf")
mygtfpath=paste0(rootpath,mygtfpath)
# output directory path for results.
#outputwd=paste0(rootpath,"bin/ggCDPbamv1results/")
outputwd=paste0("results/")
myTestName=paste0("PARCLIP_",basename(gsub(".bam","",myPARCLIPpath)),"__WT_",basename(gsub("Aligned.sortedByCoord.out.bam","",myWTpaths[1])),"__Treatment_",basename(gsub("Aligned.sortedByCoord.out.bam","",myTreatmentpaths[1])))
# PARAMETERS AND PROCESSING
myminTPM = 5

ProcessedTables <- ggCDPbamv1(gtfpath = mygtfpath,
                              myWTpaths = myWTpaths,
                              myTreatmentpaths = myTreatmentpaths,
                              myPARCLIPpath = myPARCLIPpath,
                              ReverseBamFileStrand = TRUE, # TRUE RNAseq bam files whose strand orientation needs to be reverted. FALSE: no change to strand orientation of RNAseq bam files.
                              ispairedendread = TRUE, # TRUE for paired-end bam. FALSE: single-end bam.
                              readlength = 150, # length of reads. usually 150 for paired-end and 50 for single-end
                              Only24chromosomes=TRUE, # TRUE: filter gtf for the 24 chomosomes. FALSE: include all scafolds from gtf.
                              minTPM = myminTPM, # Minimum cut-off Gene Expression (either WT or Treatment have to be higher than this)
                              maxTPM = Inf, # Maximum cut-off Expression (both WT and Treatment have to be lower than this)
                              minexon = 5, # Minimum cut-off Exon Expression (either WT or Treatment have to be higher than this)
                              minintron = 5, # Minimum cut-off Intron Expression (either WT or Treatment have to be higher than this)
                              minExtensionRatio = 0.1, # Minimum cut-off Intron Expression (either WT or Treatment have to be higher)
                              absminTPM = 1,  # bias added to both WT and Treatment Values
                              absminexon = 1,  # bias added to both WT and Treatment Values
                              absminintron = 1,  # bias added to both WT and Treatment Values
                              absminExtensionRatio = 0.01,  # bias added to both WT and Treatment Values
                              minTxSize = 500, # in nucleotides
                              minExonSize = 500, # in nucleotides
                              minIntronSize = 200, # in nucleotides
                              Extension = 1000, # length of downstream region
                              RNAseqLengthNormalize = TRUE, # TRUE: normalizing RNA-seq by length (RPKM). FALSE: no normalization (TPM).
                              TargetTolerance=50, # target proximity tolerance
                              ignoreParclipStrand = FALSE,
)

MT <- ProcessedTables[[1]]
MTreplicates <- ProcessedTables[[2]]
write.table(MT, file = paste0(outputwd,myTestName,"_MasterTable.csv"), sep = "\t", row.names = FALSE)
write.table(MTreplicates, file = paste0(outputwd,myTestName,"_MasterTablereplicates.csv"), sep = "\t", row.names = FALSE)
write.table(allresults$AveragedValues[1:10,], file = paste0(outputwd,myTestName,"_MasterTable.csv"), sep = "\t", row.names = FALSE)

allresults$AveragedValues
# check
MT[gene_id %in% "ENSG00000064932.15",]           # see results for this gene ensemble id
MT[gene_name == "SBNO2",]                        # see results for this gene name
MT[gene_name %in% c("KHSRP","FUBP1","FUBP3"),]   # see results for this list of gene names
MT[grepl("X3",gene_name),]                       # see results for all gene names that contain characters X3
MT[grepl("^FUBP",gene_name),]                    # see results for all gene names that begin with characters FUBP
MT[grepl("X3$",gene_name),]                      # see results for all gene names that end with characters X3



## PCA and Reproducibility

#PCA
# library(ggfortify)
# library(DESeq2)
genecolnames <- colnames(MTreplicates)[which(grepl("WTgene|Treatmentgene",colnames(MTreplicates)))]
sampleinfo <- data.frame(identifier=genecolnames)
sampleinfo$exp <- ifelse(grepl("WTgene",sampleinfo$identifier), 'WT', 'Treatment')

replicates <- MTreplicates[,genecolnames, with=FALSE]
#replicates[, TPM := rowMeans(replicates[,genecolnames,with=FALSE],na.rm = TRUE)]
#replicates <- replicates[TPM >= myminTPM,]
replicates <- log2(replicates + 1)
replicates.pca <- prcomp(t(as.matrix(replicates)))
replicates.pca
summary(replicates.pca)
replicates.pca.plot <- autoplot(replicates.pca,
                                data = sampleinfo, 
                                colour="exp", 
                                shape=FALSE,
                                labels="identifier",
                                size=1)
replicates.pca.plot

# HEATMAP
countVar <- apply(replicates, 1, var)
# Get the row numbers for the top 500 most variable genes
highVar <- order(countVar, decreasing=TRUE)[1:500]
# Subset logcounts matrix
hmDat <- replicates[highVar,]
mypalette <- brewer.pal(11, "RdYlBu")
# http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$exp]
# Plot the heatmap
heatmap.2(as.matrix(hmDat), 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 500 most variable genes across samples",
          ColSideColors=col.cell,scale="row",cexCol = 0.5,
          keysize=4)


myidr <- function(x,mymu,mysigma){
  mu <- mymu
  sigma <- mysigma
  rho <- 0.8
  p <- 0.5
  if (dim(x)[2]>1){
  idr.out <- est.IDR(x, mu, sigma, rho, p, eps=0.001, max.ite=30)
  # select observations exceeding IDR threshold=0.05
  #names(idr.out)
  #idr.out$para
  IDR.level <- 0.05
  x.selected <- select.IDR(x, idr.out$IDR, IDR.level)
  (x.selected$n)/(dim(x)[1])
  return((x.selected$n)/(dim(x)[1]))
  } else {return(1)
    }
}

WTgenecolnames <- colnames(MTreplicates)[which(grepl("WTgene",colnames(MTreplicates)))]
WTreplicates <- MTreplicates[,WTgenecolnames, with=FALSE]
WTreplicates[, WTTPM := rowMeans(WTreplicates[,WTgenecolnames,with=FALSE],na.rm = TRUE)]
WTreplicates <- WTreplicates[WTTPM >= myminTPM,]
WTreplicates <- log2(WTreplicates + 1)
WTIrreproducibility <- myidr(WTreplicates[,1:(dim(WTreplicates)[2]-1)],
                             mean(WTreplicates$WTTPM,na.rm = TRUE),
                             sd(WTreplicates$WTTPM,na.rm = TRUE))

Treatmentgenecolnames <- colnames(MTreplicates)[which(grepl("Treatmentgene",colnames(MTreplicates)))]
Treatmentreplicates <- MTreplicates[,Treatmentgenecolnames, with=FALSE]
Treatmentreplicates[, TreatmentTPM := rowMeans(Treatmentreplicates[,Treatmentgenecolnames,with=FALSE],na.rm = TRUE)]
Treatmentreplicates <- Treatmentreplicates[TreatmentTPM >= myminTPM,]
Treatmentreplicates <- log2(Treatmentreplicates)
TreatmentIrreproducibility <- myidr(Treatmentreplicates[,1:(dim(Treatmentreplicates)[2]-1)],
                                    mean(Treatmentreplicates$TreatmentTPM,na.rm = TRUE),
                                    sd(Treatmentreplicates$TreatmentTPM,na.rm = TRUE))
TotalIrreproducibility <- data.frame(Sample = c("WT","Treatment"), Percentage = c(WTIrreproducibility,TreatmentIrreproducibility))

xx <- barplot(TotalIrreproducibility$Percentage,
              main = "Replicate Reproducibility", names.arg = TotalIrreproducibility$Sample,
              xlab = "Samples",ylab = " Reproducibility (%), IDR level 0.05",
              col = "blue",
              ylim = c(0,1))

text(x = xx, y = TotalIrreproducibility$Percentage, label = paste(100*round(TotalIrreproducibility$Percentage,3)," %"), pos = 3, cex = 0.8, col = "black")

# SAVE RESULTS
#======================================================================
# save plots ...
plotByColumns <- c("log2fC","log2fCexon","log2fCintron")
binColumns <- c("PARCLIPgene","PARCLIP5UTR","PARCLIPCDS","PARCLIPintron","PARCLIP3UTR") 

for (myplotByColumn in plotByColumns){
  for (mybinColumn in binColumns){
    p <- ggCDPbamv1(MasterTable = MT, plotByColumn=myplotByColumn, binColumn=mybinColumn, AnalysisName = myTestName, exprNormalization = TRUE, PARCLIPLengthNormalize = TRUE, includeBin0 = TRUE)
    ggsave(paste0(outputwd,myTestName,"_", myplotByColumn,"_intersect_",mybinColumn,".pdf"), width = 15, height = 5, p, "pdf")
  }
}
# save plots  DONE
# ================

# save MasterTable BEGIN
pdf(paste0(outputwd,myTestName,"_targetgenes_histogram.pdf"))
par(mfrow=c(2,2))
allresults <- ggCDPbamv1(MasterTable = MT, plotByColumn="log2fC", binColumn="PARCLIPgene", 
                         minTPM = 5, maxTPM = Inf, absminTPM =  1, minTxSize = 500,
                         exprNormalization = TRUE, # TRUE: caliborate binning by WT Expression. FALSE: just use log2 of Normalized XL
                         PARCLIPLengthNormalize = TRUE,
                         returnAll = TRUE)
hist(allresults$AveragedValues$target, col = "blue", breaks = 100, main = "Target distribution\nwith RNAseq Calibration\n with Length Calibration")
allresults <- ggCDPbamv1(MasterTable = MT, plotByColumn="log2fC", binColumn="PARCLIPgene", 
                         minTPM = 5, maxTPM = Inf, absminTPM =  1, minTxSize = 500,
                         exprNormalization = FALSE, # TRUE: caliborate binning by WT Expression. FALSE: just use log2 of Normalized XL
                         PARCLIPLengthNormalize = TRUE,
                         returnAll = TRUE)
hist(allresults$AveragedValues$target, col = "blue", breaks = 100, main = "Target distribution\nwithout RNAseq Calibration\n with length Calibration")

allresults <- ggCDPbamv1(MasterTable = MT, plotByColumn="log2fC", binColumn="PARCLIPgene", 
                         minTPM = 5, maxTPM = Inf, absminTPM =  1, minTxSize = 500,
                         exprNormalization = TRUE, # TRUE: caliborate binning by WT Expression. FALSE: just use log2 of Normalized XL
                         PARCLIPLengthNormalize = FALSE,
                         returnAll = TRUE)
hist(allresults$AveragedValues$target, col = "blue", breaks = 100, main = "Target distribution\nwith RNAseq Calibration\n without Length Calibration")
allresults <- ggCDPbamv1(MasterTable = MT, plotByColumn="log2fC", binColumn="PARCLIPgene", 
                         minTPM = 5, maxTPM = Inf, absminTPM =  1, minTxSize = 500,
                         exprNormalization = FALSE, # TRUE: caliborate binning by WT Expression. FALSE: just use log2 of Normalized XL
                         PARCLIPLengthNormalize = FALSE,
                         returnAll = TRUE)
hist(allresults$AveragedValues$target, col = "blue", breaks = 100, main = "Target distribution\nwithout RNAseq Calibration\n without length Calibration")
dev.off()

# saving MasterTable ...
allresults <- ggCDPbamv1(MasterTable = MT, plotByColumn="log2fC", binColumn="PARCLIPgene", 
                         minTPM = 5, maxTPM = Inf, absminTPM =  1, minTxSize = 500,
                         exprNormalization = TRUE, # TRUE: caliborate binning by WT Expression. FALSE: just use log2 of Normalized XL
                         PARCLIPLengthNormalize = TRUE,
                         returnAll = TRUE)
write.table(allresults$AveragedValues, file = paste0(outputwd,myTestName,"_MasterTable.csv"), sep = "\t", row.names = FALSE)
# save MasterTable Done
#====================
# MIXTURE MODEL
#library('mclust')
#https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html
y=unlist(allresults$AveragedValues$target[which(is.finite(allresults$AveragedValues$target))])
mean(y)
yBIC = mclustBIC(y, modelNames="V")
yModel = mclustModel(y, yBIC)
yModel$parameters$mean
yModel$parameters$variance$sigmasq
hist(allresults$AveragedValues$target, col = "blue", breaks = 100, main = "Target distribution\nwithout RNAseq Calibration\n without length Calibration", freq = FALSE)
curve(dnorm(x,mean=mean(y),sd=sd(y)), add=TRUE,col="dark blue",lwd=3)
for (m in 1:length(yModel$parameters$mean)){
  curve(yModel$parameters$pro[m]*dnorm(x,mean=yModel$parameters$mean[m],sd=yModel$parameters$variance$sigmasq[m]), add=TRUE,col="red",lty=2)
}

# SEE RESULTS
#======================================================================
#See histogram of binning
allresults <- ggCDPbamv1(MasterTable = MT, plotByColumn="log2fC", binColumn="PARCLIP3UTR", 
                         minTPM = 5, maxTPM = Inf, absminTPM =  1, minTxSize = 500,
                         exprNormalization = TRUE, # TRUE: caliborate binning by WT Expression. FALSE: just use log2 of Normalized XL
                         PARCLIPLengthNormalize = TRUE,
                         includeBin0 = FALSE,
                         returnAll = TRUE)
hist(allresults$AveragedValues$target,breaks = 100,col = "blue")

# CDP plots
#Options: 
#change plotByColumns="log2fC" to plotByColumns="log2fCintron" or plotByColumns="log2fCextension")
#change binColumns="PARCLIPgene" to "PARCLIP5UTR","PARCLIPCDS","PARCLIPintron", or "PARCLIP3UTR" 

MT <- read.csv(file = paste0(outputwd,myTestName,"_MasterTable.csv"), sep="\t",header = TRUE)
MT <- setDT(MT)
ggCDPbamv1(MasterTable = MT, plotByColumn="log2fC", binColumn="PARCLIPgene",
           minTPM = 5, maxTPM = Inf, absminTPM =  1, minTxSize = 500,        # gene expression filter (always applied to any plot)
           minintron = 5, absminintron =  1, minIntronSize = 200, # intron expression filter
           minexon = 5, absminexon =  1, minExonSize = 500,
           minExtensionRatio = 0.1, absminExtensionRatio = 0.01,    # downstream expression filter
           exprNormalization = TRUE, # TRUE: calibrate binning by WT expression. False: Just use log2 of normalized XL
           PARCLIPLengthNormalize = TRUE,
           includeBin0 = TRUE, # TRUE: include bin0. FALSE: exclude bin0 and p-values.
           lineSize = 0.2,
           panelBorderSize = 0.5,
           legendSpacing = 1,
           xlab = "Log2 Fold Change [Treatment/WT]",
           ylab = "Cumulative Distribution",
           axisTitleSize = 11,
           axisTextSize = 11,
           ggTitleSize = 11,
           returnAll = TRUE
)   

# 3D plot
fig <- plot_ly(allresults$binnedTable, x = ~log2fCexon, y = ~log2fCintron, z = ~log2fC, text = ~paste0(gene_name,": ",gene_type) , color = ~binning, colors = allresults$plotColors, size = ~1, sizes = c(1,10), mode = 'text')
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'log2fCexon'),
                                   yaxis = list(title = 'log2fCintron'),
                                   zaxis = list(title = 'log2fC'))) 
fig
#use the following commented line to save 3D plot as html file
#htmlwidgets::saveWidget(as_widget(fig), paste0(outputwd,myTestName,"_3D.html"))
#MAEL
#
#myPARCLIPpath = paste0(rootpath,"Pavol/MAEL/MH55/MH55_MAEL.bam")
#myWTpaths = c("/home/manzouro/Pavol/MAEL/Hisat/HEK_NONE_DOX_rep1.bam","/home/manzouro/Pavol/MAEL/Hisat/HEK_NONE_DOX_rep2.bam","/home/manzouro/Pavol/MAEL/Hisat/HEK_NONE_DOX_rep3.bam")
#myTreatmentpaths = c("/home/manzouro/Pavol/MAEL/Hisat/HEK_MAEL_DOX_rep1.bam","/home/manzouro/Pavol/MAEL/Hisat/HEK_MAEL_DOX_rep2.bam","/home/manzouro/Pavol/MAEL/Hisat/HEK_MAEL_DOX_rep3.bam")

#KAZU
# Kazuexp236_KA_29_34_RPI12_S12_L002_R1_001.aligned_TorAwithheader_sorted.bam
# Kazuexp236_KA_29_34_RPI12_S12_L002_R1_001.aligned.bam
#myPARCLIPpath = c("/home/manzouro/Meguro/SBNO2/temp/PARCLIP/Kazuexp236_KA_29_34_RPI12_S12_L002_R1_001.aligned_TorAwithheader_sorted.bam")
#if (!RunOnServer) {
#  myPARCLIPpath=paste0(rootpath,gsub(myaccount,"",myPARCLIPpath))
#}
#myWTpaths = c("/home/manzouro/Meguro/SBNO2/final/NC_IFNb2h_Nuc1_S3_L008_R1_001Aligned.sortedByCoord.out.bam","/home/manzouro/Meguro/SBNO2/final/NC_IFNb2h_Nuc2_S6_L008_R1_001Aligned.sortedByCoord.out.bam")
#myTreatmentpaths = c("/home/manzouro/Meguro/SBNO2/final/SBNO2KD_IFNb2h_Nuc1_S21_L008_R1_001Aligned.sortedByCoord.out.bam","/home/manzouro/Meguro/SBNO2/final/SBNO2KD_IFNb2h_Nuc2_S24_L008_R1_001Aligned.sortedByCoord.out.bam")
#if (!RunOnServer) {
#  myWTpaths=paste0(rootpath,gsub(myaccount,"",myWTpaths))
#  myTreatmentpaths=paste0(rootpath,gsub(myaccount,"",myTreatmentpaths))
#}
#mygtfpath=c("/home/manzouro/data/Homo_sapiens/NCBI/GRCh38/NOALT/gencode.v25.primary_assembly.annotation.sorted.gtf")
#if (!RunOnServer) {
#  mygtfpath=paste0(rootpath,gsub(myaccount,"",mygtfpath))
#}

#Dimitrios
#myPARCLIPpath = c("/home/manzouro/temp/testFUBPKHSRP/PARCLIP/KHSRP/FUBP2_29_02_ABC.aligned_TtoC.bam")
#if (!RunOnServer) {
#  myPARCLIPpath=paste0(rootpath,gsub(myaccount,"",myPARCLIPpath))
#}


#myWTpaths = c("/home/manzouro/temp/testFUBPKHSRP/PAIREDEND/HEKAAligned.sortedByCoord.out.bam","/home/manzouro/temp/testFUBPKHSRP/PAIREDEND/HEKBAligned.sortedByCoord.out.bam","/home/manzouro/temp/testFUBPKHSRP/PAIREDEND/HEKCAligned.sortedByCoord.out.bam")
#myTreatmentpaths = c("/home/manzouro/temp/testFUBPKHSRP/PAIREDEND/FUBP2AAligned.sortedByCoord.out.bam","/home/manzouro/temp/testFUBPKHSRP/PAIREDEND/FUBP2BAligned.sortedByCoord.out.bam","/home/manzouro/temp/testFUBPKHSRP/PAIREDEND/FUBP2CAligned.sortedByCoord.out.bam")
#if (!RunOnServer) {
#  myWTpaths=paste0(rootpath,gsub(myaccount,"",myWTpaths))
#  myTreatmentpaths=paste0(rootpath,gsub(myaccount,"",myTreatmentpaths))
#}
#mygtfpath=c("/home/manzouro/data/Homo_sapiens/NCBI/GRCh38/NOALT/gencode.v25.primary_assembly.annotation.sorted.gtf")
#if (!RunOnServer) {
#  mygtfpath=paste0(rootpath,gsub(myaccount,"",mygtfpath))
#}


#myPARCLIPpath = c("FUBP2_29_02_ABC.aligned_TtoC.bam")
#if (!RunOnServer) {
#  myPARCLIPpath=paste0(rootpath,gsub(myaccount,"",myPARCLIPpath))
#}
#myWTpaths = c("HEKAAligned.sortedByCoord.out.bam","HEKBAligned.sortedByCoord.out.bam","HEKCAligned.sortedByCoord.out.bam")
#myTreatmentpaths = c("FUBP2AAligned.sortedByCoord.out.bam","FUBP2BAligned.sortedByCoord.out.bam","FUBP2CAligned.sortedByCoord.out.bam")
#mygtfpath=c("gencode.v25.primary_assembly.annotation.gtf")


