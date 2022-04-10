
# VST transformed
gtex <- readRDS("gtexv8_vst_combat3_dds-2.Rds")
test <- read.table(file='sa_location_Secreted To Blood.tsv',sep = '\t', header = TRUE)


SubsetGtexByTissue <- function(gtex, tissue)
{
  # Get a list of subject IDs with RNA-seq data in the tissue of interest
  # Args:
  #       gtex: Gtex matrix DESeq2 data set, along with sample annotation and demographics in the colData
  #       tissue: str tissue of interest corresponding to the SMTSD field of the Gtex annot, e.g., "Heart - Left Ventricle"
  # Output:
  #      List of GTEx sample
  
  coldata <- colData(gtex)
  samples <- coldata[which(coldata$SMTSD == tissue), ]
  
  return(rownames(samples))
  
}




PlotGTExData <- function(gtex, gene, tissue){
  # Plots out age-associated changes in GTEx data of a particular gene
  # 
  # Args:
  #       gtex: Gtex matrix DESeq2 data set, along with sample annotation and demographics in the colData
  #       gene: str of gene name, e.g., "KLF15
  #       tissue: str tissue of interest corresponding to the SMTSD field of the Gtex annot, e.g., "Heart - Left Ventricle"
  #
  # Output:
  #      GGplot object comparing expression across age.
  
  require(dplyr)
  require(ggplot2)
  require(rstatix)
  require(ggpubr)
  
  # Get GTEx samples that correspond to tissue
  sample_names <- SubsetGtexByTissue(gtex, tissue)
  
  # Subset Gtex matrix to gene of interest and only tissue samples
  gtex_gene <- assay(gtex)[gene, sample_names] %>% as.data.frame() %>% rownames_to_column()
  
  colnames(gtex_gene) <- c("SAMPID", "value")
  gtex_gene <- gtex_gene %>% dplyr::left_join(dplyr::select(as.data.frame(colData(gtex)), SAMPID, SUBJID, AGE))
  
  gtex_gene$AGEGROUP <- factor(gtex_gene$AGE)
  gtex_gene$AGE <- as.numeric(gsub("-.*", "", gtex_gene$AGE))
  
  # Perform ANOVA across age groups (stat.test <- aov(value ~ AGE, data = gtex_sub) %>% tukey_hsd())
  
  # GGPlot of boxplot over age
  g <-ggscatter(gtex_gene, x = "AGE", y = "value",add = "reg.line", conf.int = TRUE,add.params = list(color = "blue",
                                                                                                      fill = "lightgray"))
  g<-g+stat_cor(method = "pearson",p.accuracy = 0.01, r.accuracy = 0.01,label.x.npc = 0.05, 
                label.y.npc = 1,r.digits = 3,p.digits = 3) 
  g <- g + ggthemes::theme_few() + theme(aspect.ratio=1, legend.position="none", text=element_text(size=6), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) 
  g <- g + ylab("Normalized RNA Expression (VSD)")
  g <- g + ggtitle(gene, tissue) + ggplot2::scale_fill_viridis_d(begin=0.3, end=0.8)
  g
  
}


getGTExCorrelation <- function(gtex,test,gene,tissue)
  {
  # Get the correlation between GTEx slope and a gene list
  # 
  # Args:
  #       gtex: Gtex matrix DESeq2 data set, along with sample annotation and demographics in the colData
  #       res: A data frame with gene name from the DDS result
  #       tissue: str tissue of interest corresponding to the SMTSD field of the Gtex annot, e.g., "Heart - Left Ventricle"
  #
  # Output:
  #      A result dataframe containing log2FC in the DDS result as well as GTEx v8 correlation
  
  require(dplyr)
  require(WGCNA)
  
  
  # Get GTEx samples that correspond to tissue
  sample_names <- SubsetGtexByTissue(gtex, tissue)
  
  ddsRes_sub <- test %>% dplyr::select(Gene)
  ddsRes_sub$gtex_corr <- NA
  ddsRes_sub$gtex_pearson_p <- NA
  
  
  # Get the lsit of available GTEx genes and subset the result file with it
  gtex_available_genes <- rownames(assay(gtex))
  
  # Subset Gtex matrix to gene of interest and only tissue samples
  for(i in 1:nrow(test)){
    
    gene <- toupper(test$Gene[i])
    
    if (gene %in% gtex_available_genes){
      sample_names <- SubsetGtexByTissue(gtex, tissue)
      
      # Subset Gtex matrix to gene of interest and only tissue samples
      gtex_gene <- assay(gtex)[gene, sample_names] %>% as.data.frame() %>% rownames_to_column()
      
      colnames(gtex_gene) <- c("SAMPID", "value")
      gtex_gene <- gtex_gene %>% dplyr::left_join(dplyr::select(as.data.frame(colData(gtex)), SAMPID, SUBJID, AGE))
      
      gtex_gene$AGEGROUP <- factor(gtex_gene$AGE)
      gtex_gene$AGE <- as.numeric(gsub("-.*", "", gtex_gene$AGE))
      
      # Perform ANOVA across age groups (stat.test <- aov(value ~ AGE, data = gtex_sub) %>% tukey_hsd())
       ddsRes_sub$gtex_corr[i] <- cor(gtex_gene$AGE, gtex_gene$value,method="pearson")
       ddsRes_sub$gtex_pearson_p[i]<-cor.test(gtex_gene$AGE, gtex_gene$value,method="pearson")$p.value      
    }
    
  }
  
  return (ddsRes_sub)
}






gtex_available_genes <- rownames(assay(gtex))

tissue="Heart - Left Ventricle"
at<-getGTExCorrelation(gtex,test,gene,tissue)
at<-na.omit(at)
write.csv(at,paste(tissue, ".csv", sep=""))
myplots_Heart_Left_Ventricle<-list()
for(i in 1:nrow(test))
  {
  
  gene <- toupper(test$Gene[i])
  if (gene %in% gtex_available_genes)
    {
    message(i)
    myplots_Heart_Left_Ventricle[[i]] <- local({
      PlotGTExData (gtex, gene, tissue)
    })
  
   }
  
}
myplots_Heart_Left_Ventricle = myplots_Heart_Left_Ventricle[-which(sapply(myplots_Heart_Left_Ventricle, is.null))]

tissue="Heart - Atrial Appendage"
at<-getGTExCorrelation(gtex,test,gene,tissue)
at<-na.omit(at)
write.csv(at,paste(tissue, ".csv", sep=""))
myplots_atrial_appendage<-list()
for(i in 1:nrow(test))
{
  
  gene <- toupper(test$Gene[i])
  if (gene %in% gtex_available_genes)
  {
    message(i)
    myplots_atrial_appendage[[i]] <- local({
      PlotGTExData (gtex, gene, tissue)
    })
    
  }
  
}
myplots_atrial_appendage = myplots_atrial_appendage[-which(sapply(myplots_atrial_appendage, is.null))]


tissue="Adipose - Visceral (Omentum)"
at<-getGTExCorrelation(gtex,test,gene,tissue)
at<-na.omit(at)
write.csv(at,paste(tissue, ".csv", sep=""))
myplots_Adipose<-list()
for(i in 1:nrow(test))
{
  
  gene <- toupper(test$Gene[i])
  if (gene %in% gtex_available_genes)
  {
    message(i)
    myplots_Adipose[[i]] <- local({
      PlotGTExData (gtex, gene, tissue)
    })
    
  }
  
}
myplots_Adipose = myplots_Adipose[-which(sapply(myplots_Adipose, is.null))]




tissue="Muscle - Skeletal"
at<-getGTExCorrelation(gtex,test,gene,tissue)
at<-na.omit(at)
write.csv(at,paste(tissue, ".csv", sep=""))
myplots_Skeletal<-list()
for(i in 1:nrow(test))
{
  
  gene <- toupper(test$Gene[i])
  if (gene %in% gtex_available_genes)
  {
    message(i)
    myplots_Skeletal[[i]] <- local({
      PlotGTExData (gtex, gene, tissue)
    })
    
  }
  
}
myplots_Skeletal = myplots_Skeletal[-which(sapply(myplots_Skeletal, is.null))]


aging_gene<-test$Gene
g<-list()
for(i in 1:length(myplots_Heart_Left_Ventricle))
{
  if (myplots_Heart_Left_Ventricle[[i]]$labels$title  %in% aging_gene)
  {
    g[[i]]<-local({
      ggarrange(myplots_Heart_Left_Ventricle[[i]] ,myplots_atrial_appendage[[i]],myplots_Adipose[[i]],myplots_Skeletal[[i]], ncol = 4, nrow = 1)
       })
  }
  file_name = paste(myplots_Heart_Left_Ventricle[[i]]$labels$title, ".tiff", sep="")
  ggsave(file_name,width = 10, height = 8, dpi = 300, units = "in")
  print(g[[i]])
  dev.off()
}



