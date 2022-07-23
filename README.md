# Differentially Expressed Genes Associated with Biopsy Sites for Lung Adenocarcinoma
## Abstract:
This study analyzed molecular data samples from The Cancer Genome Atlas database for Lung Adenocarcinoma. This cancer specifically grows in mucus secreting glands of the lung. The data from the database was then analyzed utilizing DESeq analysis in order to identify differences in gene expression between upper and lower lobe lung biopsy samples in males with M0 lung adenocarcinoma. Rstudio with a series of packages including TCGAbiolinks was utilized to carry out these analyses. Additionally, a volcano plot, one-gene plot, and a principle component analysis plot was created in order to visualize the results. It was discovered that 616 or 1.96% of the genes were significant. Two of the most significant genes were SNORA73B and REG4. These genes were associated with cancer cell proliferation and often led to the advancement of cancer. However, there is no definite relationship of these genes to the upper and lower lobe sites of the lung. Furthermore, differences between these two areas could be due to confounding variables or smoking's effect on the specific tissues. It is then possible that the differences in biopsy sites for lung cancer may not be explained by the significantly differentially expressed genes that were isolated with this study. Follow up research is needed in order to identify the role that these significant genes play in cancer cells and how they can be isolated even further. This study could be beneficial to the treatment of cancers present in the investigated regions.

## Methods:
First the GDC query function within the TCGAbiolinks package was used to obtain the data sets involved for TCGA-LUAD. The parameters that were placed into this function included: project = TCGA-LUAD, data.category = "Transcriptome Profiling", data.type  = "Gene Expression Quantification", workflow.type = "STAR - Counts", experimental.strategy = "RNA-Seq", and legacy = F. These parameters helped narrow down the data that would be required for the differential expression analysis (National Cancer Institute).
	Then to prepare the data for efficient downloading, the project was divided into ten smaller groupings of barcodes. This data was then stored into countList which is a vector. Then all the elements of countList were combined to become a summarizedExperiment object called luad.counts.
	After all the data for TCGA-LUAD has been obtained, the data was then subsetted based on the controls which were gender=”male” and ajcc_pathologic_m=”M0”. Then all the values that were not “Lower lobe, lung” and “Upper lobe, lung” were removed from the site_of_resection_or_biopsy column, since we are only comparing lower lobe and upper lobe of the lung biopsy samples.
	The luadSubset was then subsetted even further by the mean and standard deviation of the count. A column was created for the mean and standard deviation in the subset. Then in order to visualize the mean counts and standard deviation, a scatter plot was created. Then in order to remove any extreme outliers in the data, the rows with a standard deviation of greater than 13,000 were removed from the subset. Additionally, the genes that were greater than fifty percent of the samples at zero counts were removed from luadSubset. Lastly, a column for comparison was created as “comp”, This allowed the code to successfully compare and contrast the gene expression between the upper lobe and lower lobe lung biopsy sites. The “Upper lobe, lung” was True or the manipulated variable and the “Lower lobe, lung” biopsy site was false or the control variable.
	Next, a DESeq analysis was conducted in order to compare differential gene expression between the upper lobe and lower lobe biopsy sites. The function, DESeq, was able to calculate the fold changes and p-values that help determine the significance of the gene expressions. The design parameter in the function was equal to the comp column that was created earlier in the process and the other parameters in the function were defaulted. Then the results with the gene data were placed into a vector called resOutput (Anders, 2021). 
	Finally, in order to visualize the results of the analysis, a volcano plot was created. The cutoffs for the graph were set at 0.01 and 1.50. These are the maximum false discovery rate p-value and the minimum log fold difference (cite). Then the genes were identified as either above or below the cutoffs. The volcano plot used the ggplot2 function in order to plot the different genes according to their log2 fold change and the adjusted p-values. Then the points on the volcano plot were differentiated based on the biopsy site of the patient and the level of gene expression. The top ten genes were labeled, which were the most significant in differential expression. Furthermore, one gene was chosen to compare its distribution across both groups using ggplot, as given in figure 2 (Ggplot2, 2021). Lastly, a PCA plot was created through further analysis on the data in order to visualize which samples were extremely different (Nolan Bentley, personal communication).



## Results:
In the DESeq analysis, the upper lobe biopsy samples (manipulated) and the lower lobe biopsy samples (control) were compared. This analysis established that 317 genes were statistically significant and downregulated, while 299 genes were significantly upregulated as compared to the other genes. The volcano plot as given in figure 1 can be interpreted by finding genes that are the most upregulated to be located on the right side of the plot and the most downregulated genes will be toward the left side of the plot. Additionally, the genes that are significant will be represented by a triangle as given in the Sig. DE legend in figure 1.

###### Figure 1: The volcano plot below displays the differential gene expression of genes in upper lobe and lower lobe lung cancer samples. The significant DE data is represented by a triangle and the non-significant genes are represented by a circle.  
![Volcano Plot](https://github.com/AreesaM/CancerProject/blob/aee6b32af73f6d6685c3c2aa1cfa0c24a48e1b64/Figure_Volcano1.png?raw=true)



###### Figure 2: The one-gene plot below displays the single gene distribution for the REG4 gene for both the upper lobe(TRUE)  and lower lobe (FALSE) biopsy site samples. This gene was found to be the most down-regulated gene for the manipulated group or in our case the upper lobe biopsy samples. 
![REG4](https://github.com/AreesaM/CancerProject/blob/aee6b32af73f6d6685c3c2aa1cfa0c24a48e1b64/Figure_REG4.png?raw=true)
###### Figure 3: The principle component analysis plot below displays which samples should be removed from the study, as the samples to the very right are extremely different from the localized area of the other samples in the plot. The red dots represent the lower lobe biopsy site samples and the blue dots represent the upper lobe biopsy samples. 
![PCA Plot](https://github.com/AreesaM/CancerProject/blob/4b99a50019cd86d728de84db4a22b20791867b95/pca_plot.jpg?raw=true)
As given in table 2 below, we determined that 30,753 genes or 98.04 % of the genes analyzed were not significant and 616 or 1.96% of the genes were significant. This data was determined based on if the gene fell above or below the p-value and log fold cutoffs, which represent the significance and regulation of the gene. Additionally, according to table 2, the first “FALSE” means the gene is not significantly expressed and if there is a “TRUE” it is significantly expressed. The second “FALSE” means the gene is downregulated and if it is a “TRUE” the gene is upregulated. Overall, 30,753 or 98.04% of the genes were insignificant, and 616 or 1.96% of the genes were significantly differentially expressed. 




###### Table 1: This table provides a count of genes and samples before and after subsetting.

|               |Before Subsetting|After Subsetting|
| ------------- |:---------------:|:---------------:|
|Number of Genes| 60660     |   31369 |
|Number of Lower Lobe Biopsy samples| 196     |  62 |
|Number of Upper Lobe Biopsy samples|357|101



###### Table 2: This table provides the number of genes that were either significantly downregulated or upregulated

|pAdj_log2Change| Number of Genes|
|:--------------:|:---------------:|
|FALSE FALSE|30641
|FALSE TRUE|112|
TRUE FALSE|317
TRUE TRUE|299





## Analysis:
Based on table 2, 317 genes were significant and downregulated and 299 genes were significantly upregulated. Based on these numbers, <2% of genes were differentially expressed for the upper lobe and lower lobe biopsy samples for men with M0 cancer. This is mostly likely due to the fact that only men who had M0 cancer or cancer that was localized in only one area of the body were used to complete this analysis. However, differences in gene expression for the upper lobe and lower lobe could be attributed to the differences in tissue and the overall differences of how smoking affects these areas of the lungs. There could also have been many confounding variables that could have created differences in gene expression such as stage of cancer or even age that we did not take into account. 
	In order to understand how these genes play a role in lung adenocarcinoma in the upper lobe and lower lobe regions, we can look at the most significant genes that can be seen in the volcano plot. In this analysis REG4 was the most significantly downregulated gene for the upper lobe biopsy site. According to a study the REG4 gene plays a significant role in inhibiting apoptosis of cancer cells and promotes cancer cell proliferation (Zhang, 2021). Based on this, we can assume that the upregulation of REG4 is associated with advancing tumor growth for lung adenocarcinoma. In comparison, a gene that was significantly upregulated for the manipulated or upper lobe biopsy site was the SNORA73B gene. According to a study, if one is deficient in SNORA73, it will lead to metabolic reprogramming of cells. In our case this will most likely lead to an increase in lung cancer cell proliferation as the cell cycle will be disrupted (Sletten.2021). However, there is no study that links any significant relationship between these two genes to the upper and lower lobe lung sites for lung cancer. Further research that isolates the relationship between the location of the biopsy site and the significant genes must be conducted in order to determine why these genes are differentially expressed. Therefore, we can conclude that the differences in biopsy sites for lung cancer may not be explained by the significantly differentially expressed genes that were isolated with this study.














## Works Cited

Anders, Simon, and Michael Love. Analyzing RNA-Seq Data with DESeq2. 26 Oct. 2021, 
: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential- expression-analysis.
    
Bentley, Nolan. Personal Communication. Apr. 2022.

“Encyclopedia.” National Cancer Institute, https://docs.gdc.cancer.gov/Encyclopedia/pages/HTSeq-Counts/. Ggplot2. (n.d.). Retrieved April 20, 2022, from https://ggplot2.tidyverse.org/

Sletten, A.C., Davidson, J.W., Yagabasan, B. et al. Loss of SNORA73 reprograms cellular metabolism and protects against steatohepatitis. Nat Commun 12, 5214 (2021). https://doi.org/10.1038/s41467-021-25457-y

Slowikowski, K. (n.d.). Ggrepel. Retrieved April 20, 2022, from https://ggrepel.slowkow.com.

“The Cancer Genome Atlas Program.” National Cancer Institute, https://www.cancer.gov/about-
: nci/organization/ccg/research/structural-genomics/tcga.

Zhang, J., Zhu, Z., Miao, Z., Huang, X., Sun, Z., Xu, H., & Wang, Z. (2021). The Clinical Significance and Mechanisms of REG4 in Human Cancers. Frontiers in oncology, 10, 559230. https://doi.org/10.3389/fonc.2020.559230


