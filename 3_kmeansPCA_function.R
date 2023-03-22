########################code for determining k-means groups########################
##Barbara Langille
#March 2023

##Use only the outlier region identified by lostruct 
#will need a ped/map and fam files




#load these packages
library(pcadapt)
library(gplots)
library(ggplot2)
library(readr)
library(ggman)


#get into the working directory
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())

#PCA function - shouldn't need to change anything except maybe value of K
PCA_getgroups <- function(PEDfile, FAMfile, MAPfile) {
  #get the data and metadata
  PED<- read.pcadapt(PEDfile, type = "ped")
  FAM <- read.delim(FAMfile, stringsAsFactors = F, header = F, sep = "")
  MAP <- read.delim(MAPfile, stringsAsFactors = F, header = F)
  colnames(MAP) <- c("Chrom", "SNP", "CM", "BP")
  #get the PCS for 10 components - usually more than enough
  PCs <- pcadapt(PED, K = 10, min.maf = 0.01)
  PVALS <- PCs$pvalues
  PCMAP <- as.data.frame(cbind(MAP, PVALS))
  ##calculate EV-proportion of explained variance so that I can add it to the PCA
  EV <- ((PCs$singular.values^2)*100)
  EV_sig <- formatC( signif(EV, digits=5), big.mark=",", format="fg")
  EV1 <- EV_sig[[1]][1]
  EV2 <- EV_sig[[2]][1]
  EV1[] <- paste0(as.matrix(EV1), '%')
  EV1[] <- paste0(as.matrix(EV1), '  PC1')
  EV2[] <- paste0(as.matrix(EV2), '%')
  EV2[] <- paste0(as.matrix(EV2), '  PC2')
  PC_FAM <- cbind(PCs$scores, FAM)
  write_tsv(PC_FAM, "PCscores_K10.tsv")
  #make the plot and add in what your explained variance is for each axis
  PCA_plot <- plot(PCs, option="scores", pop =  FAM$V1, i = 1, j =2) +#Change to Metadata$Pop to get colored points by population
    xlab(EV1) +
    ylab(EV2) + 
    theme_classic()
  #Plot loadings
  ggman_plot <- ggman(PCMAP, chrom = "Chrom", 
                      pvalue = "PVALS", 
                      snp = "SNP", 
                      bp="BP", 
                      pointSize = 1, 
                      title = "Loadings of PCA pvalues", 
                      xlabel = "Chromosome", 
                      ymax = 100 ) + 
    theme_classic()
  ggsave("PCA_plot.pdf", PCA_plot, width = 8, height = 6)
  ggsave("loadings_plot.pdf", ggman_plot, width = 10, height = 4)
  
}


#this is the place to put your files -PEDfile, and MAPfile
PCA_getgroups(PEDfile, FAMfile, MAPfile)


#If your plot has been working in the past and is now throwing the invalid graphics 
#state error, try resetting the graphics device by calling dev.off()
dev.off()

################################################################################





####K-means clustering based on PCA scores####




#load the following packages
library(tidyverse)
library(ggplot2)
library(readr)
library(ggforce)



#get into the working directory
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())


#read in the tsv data generated in the code above
pc_scores <- read_tsv("PCscores_K10.tsv")
#rename columns
colnames(pc_scores)[1] <- "PC1" 
colnames(pc_scores)[2] <- "PC2"
colnames(pc_scores)[3] <- "PC3"
colnames(pc_scores)[4] <- "PC4"
colnames(pc_scores)[5] <- "PC5"
colnames(pc_scores)[6] <- "PC6"
colnames(pc_scores)[7] <- "PC7"
colnames(pc_scores)[8] <- "PC8"
colnames(pc_scores)[9] <- "PC9"
colnames(pc_scores)[10] <- "PC10"
colnames(pc_scores)[11] <- "FamID"
colnames(pc_scores)[12] <- "SampleID"

#look at the PCS and pick the plot with the best separation - those PCs are the
#ones you will use moving forward
plot(pc_scores[,1:6]) 

#work off of only PC1 and 2 - could keep as many columns as you like
keep = c("PC1", "PC2")
pc_only = pc_scores[keep]

#count the within group sum of squares for different values of k 
#using from 1-20 - which is usually overkill
within_sum_squares = c()
for (i in 1:20){
  within_sum_squares[i] = sum(kmeans(pc_only,
                                     centers=i, nstart=25, 
                                     iter.max=1000)$withinss)
}

#get a plot showing the best K values
plot(1:20, within_sum_squares, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares", 
     main = "k-means clustering - Determining optimal k")

#choose the K based on the skree plot (inflection point)
optimal_k = kmeans(pc_only, 3, nstart=25, iter.max=1000)

# assign the clusters back to the dataframe, note the optimal_k object has other
#useful information too,such as the cluster means and the number of individuals 
#per cluster.
pc_scores$cluster_k3 = as.factor(optimal_k$cluster)

#plot the PCs with the cluster info
pcs_clusters <- ggplot() + geom_point(data = pc_scores, 
                                     aes(x = PC1, y = PC2, colour=cluster_k3))+
  geom_mark_ellipse(data=pc_scores, aes(x = PC1, y = PC2, fill = cluster_k3,
                                       label = cluster_k3))+
  labs(title =  "k-means clustering of PCs, 3 groups")+
  xlab("PC1 (%)") +
  ylab("PC2 (%)") ##can get these numbers from the PCA section above

#save the plot as an image
ggsave("kmeans.pdf", pcs_clusters, width = 8, height = 6)




