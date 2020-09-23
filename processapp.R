
library(shiny)
options(shiny.maxRequestSize = 3000*1024^2)

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
#  suppressWarnings(suppressMessages(library(idr)))
  suppressWarnings(suppressMessages(library(htmlwidgets)))
#  suppressWarnings(suppressMessages(library(ggfortify)))
#  suppressWarnings(suppressMessages(library(mclust)))
  
source(paste0("ggCDPbamv1.R"))

library(shinythemes)

ui <- fluidPage(theme = shinytheme("cerulean"),
  titlePanel("Input Paramters"),
  
  fluidRow(
    column(3,
           fileInput("myWTpaths", h6("WT bam files"), multiple = TRUE)),
    column(3,
           fileInput("myTreatmentpaths", h6("Treatment bam files"), multiple = TRUE)),
     column(3,
           fileInput("myPARCLIPpath", h6("PAR-CLIP bam file"))),
    column(3,
           fileInput("gtfpath", h6("GTF file")))
  ),
  
  fluidRow(
    column(2,
           h6("Reverse RNAseq orientation"),
           checkboxInput("ReverseBamFileStrand", "Reverse", value = TRUE)),
    column(2,
           h6("Paired-End"),
           checkboxInput("ispairedendread", "Paired-end", value = FALSE)),
    column(2, 
           numericInput("readlength", 
                        h6("Approx. RNAseq read Length"), 
                        value = 50)),
    column(2, 
           numericInput("ThreeutrExtension", 
                        h6("Length of 3utr extention"), 
                        value = 1000)),
    
       column(2,
           h6("TPM/RPKM"),
           checkboxInput("RNAseqLengthNormalize", "RPKM", value = FALSE))
    ),
    
  fluidRow(
    column(3, 
           numericInput("TargetTolerance", 
                        h6("PAR-CLIP gene target tolerance"),value = 50)), 
       column(3,
           h6("Filter Chromosomes"),
           checkboxInput("Only24chromosomes", "Only 24 Chromosomes", value = TRUE)),
 column(2,
           h6("Run"),
           actionButton("submit","Submit"))
  ),
    
   fluidRow(
      textOutput("MasterTableName"),
      textOutput("MasterTableReplicatesName")
   ),
 
 fluidRow(
   column(12,
          tableOutput("MT"))
 )
  
)

# Define server logic ----
server <- function(input, output) {

  observeEvent(input$do, { 
  output$MT <-  renderTable({
  inFile <- input$gtfpath$datapath
    if (is.null(inFile))
      return("NA")

#myPARCLIPpath = paste0(rootpath,"HNRNPK_P_Cyt_224_226_RPI36.aligned_TtoC.sorted.bam")
# myPARCLIPpath=paste0(rootpath,gsub(myaccount,"",myPARCLIPpath))

#myWTpaths = c("AF_P_mCh_V_CKDL200153406-1a-DY0088-AK1682_HC3L5BBXX_L7Aligned.sortedByCoord.out.bam")
#myWTpaths = paste0(rootpath,myWTpaths)
#myTreatmentpaths = c("AF_KO1_mCh_V_CKDL200153406-1a-DY0088-AK1544_HC3L5BBXX_L7Aligned.sortedByCoord.out.bam")
#myTreatmentpaths = paste0(rootpath,myTreatmentpaths)
#mygtfpath=c("genes_no_mir_snor_hist.gtf")
#mygtfpath=paste0(rootpath,mygtfpath)
# output directory path for results.
#outputwd=paste0(rootpath,"bin/ggCDPbamv1results/")
#outputwd=paste0(rootpath,"ggCDPbamv1results/")
#myTestName=paste0("PARCLIP_",basename(gsub(".bam","",myPARCLIPpath)),"__WT_",basename(gsub("Aligned.sortedByCoord.out.bam","",myWTpaths[1])),"__Treatment_",basename(gsub("Aligned.sortedByCoord.out.bam","",myTreatmentpaths[1])))
# PARAMETERS AND PROCESSING

ProcessedTables <- ggCDPbamv1(gtfpath = input$gtfpath$datapath,
                 myWTpaths = input$myWTpaths$datapath,
                 myTreatmentpaths = input$myTreatmentpaths$datapath,
                 myPARCLIPpath = input$myPARCLIPpath$datapath,
                 ReverseBamFileStrand = input$ReverseBamFileStrand, # TRUE RNAseq bam files whose strand orientation needs to be reverted. FALSE: no change to strand orientation of RNAseq bam files.
                 ispairedendread = input$ispairedendread, # TRUE for paired-end bam. FALSE: single-end bam.
                 readlength = input$readlength, # length of reads. usually 150 for paired-end and 50 for single-end
                 Only24chromosomes=input$Only24chromosomes, # TRUE: filter gtf for the 24 chomosomes. FALSE: include all scafolds from gtf.
                 minTPM = 0.01, 
                 maxTPM = Inf, 
                 minexon = 5, 
                 minintron = 5, 
                 minExtensionRatio = 0.1, 
                 absminTPM = 0.01,  
                 absminexon = 1,  
                 absminintron = 1,  
                 absminExtensionRatio = 0.05,  
                 minTxSize = 500, 
                 minExonSize = 500, 
                 minIntronSize = 200, 
                 Extension = input$ThreeutrExtension, 
                 RNAseqLengthNormalize = input$RNAseqLengthNormalize, 
                 TargetTolerance=input$TargetTolerance,
                 ignoreParclipStrand = FALSE
)

MT <- ProcessedTables[[1]]
MTreplicates <- ProcessedTables[[2]]
MT[1:20,]
## Save Processing
#write.table(MT, file = paste0(outputwd,myTestName,"_MasterTable.csv"), sep = "\t", row.names = FALSE)
#write.table(MTreplicates, file = paste0(outputwd,myTestName,"_MasterTablereplicates.csv"), sep = "\t", row.names = FALSE)
})

  output$MasterTableName <- renderText({
    myTestName=paste0("PARCLIP_",basename(gsub(".bam","",input$myPARCLIPpath$datapath)),"__WT_",basename(gsub("Aligned.sortedByCoord.out.bam","",input$myWTpaths$datapath[1])),"__Treatment_",basename(gsub("Aligned.sortedByCoord.out.bam","",input$myTreatmentpaths$datapath[1])))
    paste0("Output MasterTable name is:", myTestName,"_MasterTable.csv")
  })
  output$MasterTableReplicatesName <- renderText({
    myTestName=paste0("PARCLIP_",basename(gsub(".bam","",input$myPARCLIPpath$datapath)),"__WT_",basename(gsub("Aligned.sortedByCoord.out.bam","",input$myWTpaths$datapath[1])),"__Treatment_",basename(gsub("Aligned.sortedByCoord.out.bam","",input$myTreatmentpaths$datapath[1])))
    paste0("Output MasterTableReplicates name is:", myTestName,"_MasterTablereplicates.csv")
  })
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)


