# CDP

(to test the R shiny web app go to bottom)
To explore the R shiny seb server for CDP, download the sample Table (myMasterTable.csv) and uploaded it here (https://amanzour.shinyapps.io/mytestshiny/)

Complete Explanation:
The program intersects RNAseq (WT & KD) and PAR-CLIP data and creates Cummulative Distribution Plots for mRNA differential expression. Genes are categorized to targets and non-targets. The degree to which each gene might be a target depends on abundance of cross-linked reads on its mRNA (PAR-CLIP) in the context of its expression (WT RNA-seq). 
![picture](CDP.png)

Instructions: 
Construct Table
All you need is three files (ggCDPbamv1.R, processapp.R, and app.R). If you have Rstudio, create a shiny app project (ex. myshiny1) and copy the ggCDPbamv1.R and processapp.R files into it. Then press "run app" after opening processapp.R (you can rename processapp as app.R in the project, to have consistent naming). Input RNAseq bam files (when uploading, you can press and hold control botton to import multiple replicates at the same time), a PAR-CLIP or an e-CLIP bam file, and gene annotation preferrably gencode format. If your files are large, you can alter the processapp.R file by alowaying more memory by changing line : options(shiny.maxRequestSize = 3000*1024^2). Save Table after running.

Generate Plots:
Create another shiny app project (myshiny2) and copy ggCDPbamv1.R and app.R files in it. hen press "run app" after opening app.R. You can upload the table that was made from myshiny1.

To explore the R shiny web app (myshiny2), download the sample Table (myMasterTable.csv) and uploaded it here (https://amanzour.shinyapps.io/mytestshiny/)
Email me at amir.manzour@gmail.com if you have any questions.
