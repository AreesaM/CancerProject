projectDir <- "~/Bio321G/project2/"
Dir <- paste0(projectDir,"Resources/")
dir.create(Dir,recursive = T)
setwd(projectDir)


Project <- paste0(Dir, "Project.Rdata")
#This code searches for TCGA-LUAD in The Cancer Genome Atlas based on the data.category, data.type, workflow.type, experimental strategy, and legacy.
if (!file.exists(Project)) {
  # Sometimes you get the message 
  # "GDC server down, try to use this package later",
  # but seconds later the server works. Just keep trying.
  project <- TCGAbiolinks::GDCquery(
    project       = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    experimental.strategy = "RNA-Seq",
    legacy        = F
  )
  
  save(project, file = Project)
}else{
  load(Project, verbose = T)
  
}

#The barcode represents an individual case in TCGA-LUAD. This is a way to split and prepare the data to be downloaded.
barcode <- project$results[[1]]$cases

#I split it into 10 chunks of barcode. 
nSplit <- 10 
bcSplitStarts <- round(seq(1, length(barcode) / nSplit * (nSplit - 1), length.out = nSplit))
bcSplitEnds   <- c(bcSplitStarts[2:length(bcSplitStarts)] - 1, length(barcode))
print(cbind(bcSplitStarts,bcSplitEnds))

savedListName <- paste0(Dir, "qncList_len", nSplit, ".Rdata")

#Checks to see if the file exists, if not it will run the code.
if (!file.exists(savedListName)) {
  #Create empty lists to store split objects
  queryList <- list()
  countList <- list()
  
  #Loop over the split barcode groups
  for (i in 1:length(bcSplitStarts)) {
    tictoc::tic() #To allow for tracking of script timing
    #Extract current range
    currBcRange     <- bcSplitStarts[i]:bcSplitEnds[i]
    currentBarcodes <- barcode[currBcRange]
    
    #Search for the current query files
    currQuery <- GDCquery(
      project       = "TCGA-LUAD",
      data.category = "Transcriptome Profiling",
      data.type     = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      experimental.strategy = "RNA-Seq",
      legacy        = F,
      barcode       = currentBarcodes
    )
    
    #Download the current query files
    GDCdownload(currQuery, directory = paste0(projectDir, "GDCdata/"))
    
    #Prepare the current query files
    currentCount <- GDCprepare(currQuery, directory = paste0(projectDir, "GDCdata/"))
    
    #Save the output into list objects
    queryList[[i]] <- currQuery
    countList[[i]] <- currentCount
    
    #Print progress
    print(paste0("nSplit ", i, " done!"))
    tictoc::toc()
  }
  #Save the list files
  save(queryList,countList,file = savedListName)
}else{
  #Load the list object
  load(savedListName,verbose = T)
}
countList
#This combines the lists within countList into an object.
luad.counts <- do.call(cbind,countList)
#Saves file to my computer.
save(luad.counts, file = paste0(Dir,"TCGA-LUAD_RNASeq_rawCounts.rda"))
rm(list = c("queryList","countList"))
print(luad.counts)

luadDataFile <- paste0(Dir, "TCGA-LUAD_RNASeq_rawCounts.rda")
load(file = luadDataFile, verbose = T)
class(luad.counts)
mode(luad.counts)
summary(luad.counts)
print(luad.counts)
sum(table(rowRanges(luad.counts)))

##------SUBSET based on "male" and ajcc_pathologic_m status
luadSubset <- luad.counts[, which(luad.counts$gender == "male")]
luadSubset <- luadSubset[, !is.na(luadSubset$ajcc_pathologic_m)]
luadSubset <- luadSubset[,luadSubset$ajcc_pathologic_m=="M0"]
biopsylogic <-which(luadSubset$site_of_resection_or_biopsy=="Lower lobe, lung" | luadSubset$site_of_resection_or_biopsy=="Upper lobe, lung")
luadSubset<-luadSubset[,biopsylogic]
luadDf <- View(as.data.frame(colData(luadSubset)))

#This counts the number of biopsy samples before and after subsetting.
table(luad.counts$site_of_resection_or_biopsy)
table(luadSubset$site_of_resection_or_biopsy)


##---------------ANALYSIS
## Subset based on sequencing statistics
#Counts the mean counts in each of the individual cases from each gene.
luadSubset$meanCnt <- colMeans(assays(luadSubset)[[1]], na.rm = T)
#Finds the standard deviation of the counts of the cases.
luadSubset$sdCnt   <- apply(assays(luadSubset)[[1]], 2, sd, na.rm = T)
#This plots a scatter plot of the mean counts and standard deviation. It plots a red line at 13000 to remove any very different cases.
plot(luadSubset$meanCnt, luadSubset$sdCnt); abline(h = 13000,col = "red")

#Removes the cases that had a standard deviation over 13000
luadSubset <- luadSubset[, luadSubset$sdCnt < 13000]

# Subset the genes to those with less than or equal to 50% of the samples at 0 counts
## Find the mean frequency of counts equal to 0
isCnt0   <- assays(luadSubset)[[1]] == 0
meanCnt0 <- rowMeans(isCnt0)
hist(meanCnt0,1000); abline(v = 0.5,col = "red")
#I found the mean frequency of counts equal to 0
#Then I subset based on the mean frequency less than or equal to 50%

## Subset based on the mean frequency
luadSubset <- luadSubset[meanCnt0 <= 0.5, ]

# Create a column for comparison
## Create column
table(luadSubset$site_of_resection_or_biopsy,useNA = "always")
luadSubset$comp <- luadSubset$site_of_resection_or_biopsy == "Upper lobe, lung"

# Check columns
table(luadSubset$comp)
class(luadSubset$comp)

# Make a factor
luadSubset$comp <- as.factor(luadSubset$comp)


# Run DESeq
## Calculate the fold changes and p-values
ddsSE <- DESeqDataSet(luadSubset, design = ~ comp) # Converts to DESeq input
dds   <- DESeq(ddsSE) # Calculates the DE analysis
res   <- results(dds) # Organizes the results


## Add in gene data
resOutput <- cbind(as.data.frame(rowRanges(luadSubset)),as.data.frame(res))

# Set and check cutoffs
## Maximum false discovery rate adjusted p-value cutoff:
fdr.cut.off <- 0.01
## Minimum Log fold DIFFERENCE (absolute value of change) cutoff to be considered 
lfc.cut.off <- 1.50

## Check cutoffs
resOutput$padjLteCutoff    <- resOutput$padj                <= fdr.cut.off
resOutput$absFoldGteCutoff <- abs(resOutput$log2FoldChange) >= lfc.cut.off

resOutput$sigDE <- resOutput$absFoldGteCutoff & resOutput$padjLteCutoff

## Summarize results
table(paste0("pAdj_log2Change:",resOutput$padjLteCutoff,"_",resOutput$absFoldGteCutoff))
table(resOutput$sigDE)
table(resOutput$sigDE)/sum(table(resOutput$sigDE))

# Visualizing the volcano plot 
resPlotted <- resOutput[sample(1:nrow(resOutput),size = nrow(resOutput)),] # Randomize order to avoid positional biases
volc1 <- ggplot(resPlotted, aes(x     = log2FoldChange,
                                y     = -log10(padj),
                                color = baseMean ,
                                shape = sigDE)) +
  geom_point() +
  geom_vline(aes(xintercept = -lfc.cut.off)       , color = "blue" , linetype = "dashed") +
  geom_vline(aes(xintercept =  lfc.cut.off)       , color = "red"  , linetype = "dashed") +
  geom_hline(aes(yintercept = -log10(fdr.cut.off)), color = "black", linetype = "dashed") +
  labs(title = "TCGA-LUAD, Differential Gene Expression\nBetween Upperlobe and lowerlobe lung adenocarcinoma ",
                  
       y     = "-log10(FDR-Adjusted p-values)",
       x     = "log2(Fold Change)",
       color = "Base mean count",
       shape = "Sig. DE") +
  theme_bw() +
  geom_text_repel(data = head(resPlotted[order(resPlotted$padj),],10),
                  aes(label = gene_name),
                  show.legend = FALSE,color = "black",alpha = 0.8) +
  scale_color_viridis_c(trans ="log10")

volc1


##
#Investigate single row
geneInv.index    <- which.min(resPlotted$padj)
geneInv.name     <- resPlotted$gene_name[geneInv.index]
geneInv.ensemble <- rownames(resPlotted)[geneInv.index]
normCounts <- as.data.frame(counts(dds, normalized = T))
geneInv<- data.frame(
  sample     = rownames(colData(dds)),
  group      = luadSubset$comp,
  normCounts = unlist(normCounts[match(geneInv.ensemble,rownames(normCounts)),])
)

#Visualize a single gene based on grouping
##Rank the counts based on grouping
geneInv$rank <- rank(geneInv$normCounts,ties.method = "random")
##Scale from 0 to 1
geneInv$rank <- (geneInv$rank-min(geneInv$rank)) / (max(geneInv$rank) - min(geneInv$rank))
##Scale from -0.4 to 0.4
geneInv$rank <- geneInv$rank*0.8 - 0.4
##Add integer based on grouping to stagger
geneInv$rankWithFactor <- geneInv$rank + as.numeric(geneInv$group)

#Visualize single gene distribution
oneGenePlot <- ggplot(geneInv, aes(group, log10(normCounts+1)))+
  geom_boxplot(outlier.shape = NA,color = "grey50")+
  geom_point(mapping = aes(x = rankWithFactor),
             size = 3,shape = 8)+
  labs(title=paste0(geneInv.name))+
  scale_color_viridis_c()+theme_bw()

oneGenePlot
#PCA of samples
pca <- prcomp(t(normCounts),scale. = T)
plot(pca$x[,1],pca$x[,2])

# Re-plot this using an interactive plot
## This reveals that some additional samples should have been removed becuase they are extremely different!
plottedDf<-data.frame(pca$x)
plottedDf <- cbind(plottedDf,colData(dds))
scatterD3::scatterD3(PC1,PC2,data = plottedDf,
                     lab = sample_submitter_id,
                     labels_size = 0,
                     col_var = comp,
                     fixed = T)

# Output summaries
ggsave(filename = "Figure_Volcano1.png",volc1,
       dpi = 600, width = 84*2,height = 84,units = "mm")
ggsave(filename = paste0("Figure_",geneInv.name,".png"),oneGenePlot,
       dpi = 600, width = 84*2,height = 84,units = "mm")
write.csv(table(resOutput$padjLteCutoff,resOutput$absFoldGteCutoff),
          file = "Table_DEAnalysisSummary.csv")
write.csv(resOutput[order(resOutput$padj,resOutput$gene_name),],
          file = "SupplementalTable_DEAnalysisResults.csv")
write.csv(colData(luadSubset),file = "SupplementalTable_SampleMetadata.csv")






