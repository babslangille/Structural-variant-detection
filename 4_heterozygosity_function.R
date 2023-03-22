######################        heterozygosity     #############################
##Barbara Langille March 2023

##this function assesses heterozygosity of the outlier region
#ideally you have distinct clusters (ideally 3) in which you can assign individuals back to (from k-means)




#load these packages
library(dplyr)
library(ggplot2)


#Set the working directory for this script
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())



#Heterozygosity by genotype in plink - need to have a keep file of the individuals that belong to each group
#plink --file FILENAME  --hardy --maf 0.01 --keep Hap1_ch2.txt --recode  --out CHR2_Hap1 --chr 2
#plink --file FILENAME  --hardy --maf 0.01 --keep Hap2_ch2.txt --recode  --out CHR2_Hap2 --chr 2
#plink --file FILENAME  --hardy --maf 0.01 --keep Het_ch2.txt --recode  --out CHR2_Het --chr 2


#Take files from plink and load them in
Hap1  <- read.csv("CHR2_Hap1.hwe", sep="")
Hap2  <- read.csv("CHR2_Hap2.hwe", sep="")
Het  <- read.csv("CHR2_Het.hwe", sep="")

#Get the map file
MAP  <- read.csv("CHR2.map", sep="", header = F)
colnames(MAP) <- c("Chrom", "SNP", "CM", "BP")
Hap1 <- inner_join(Hap1, MAP)
Hap2 <- inner_join(Hap2, MAP)
Het <- inner_join(Het, MAP)

#make a plot with all three groups colored differently
ggplot() + geom_smooth(data = Hap1, aes(x = BP, y = O.HET., colour = "HAP1"), 
                       method = "loess",span = 0.03, se = F ) +
  geom_smooth(data = Hap2, aes(x = BP, y = O.HET., colour = "HAP2"), 
              method = "loess",span = 0.03, se = F ) +
  geom_smooth(data = Het, aes(x = BP, y = O.HET., colour = "HET"), 
              method = "loess",span = 0.03, se = F ) + geom_vline(xintercept = 5750000) + 
  geom_vline(xintercept = 12000000) + theme_classic()




################################################################################
################################################################################



#make a box plot with just the outlier region


#get a file of the observed heterozygosity and clusters (Hap1=1, Het=2, Hap2=3)
#Run the plink code above, except dont use chr function, 
#instead use a file with the snps to extract from outlier region
#Make sure the hwe file has only one tab between columns

het_function <- function(HAP1hwe, HAP2hwe, HEThwe){
  #Take files from plink and load them in
  Hap1_outliers  <- read.csv(HAP1hwe, sep="")
  Hap2_outliers  <- read.csv(HAP2hwe, sep="")
  Het_outliers  <- read.csv(HEThwe, sep="")
  #Add new col with the cluster number
  Hap1_newCol <- Hap1_outliers %>% mutate(newCol = 1)
  Hap2_newCol <- Hap2_outliers %>% mutate(newCol = 3)
  Het_newCol <- Het_outliers %>% mutate(newCol = 2)
  #get just the 2 rows we want (O.HET and newCols)
  Hap1_reduced <- Hap1_newCol %>% select(-1,-2,-3,-4,-5,-6,-8,-9)
  Hap2_reduced <- Hap2_newCol %>% select(-1,-2,-3,-4,-5,-6,-8,-9)
  Het_reduced <- Het_newCol %>% select(-1,-2,-3,-4,-5,-6,-8,-9)
  #put all three files together into 2 columns
  ohet <- rbind(Hap1_reduced, Het_reduced, Hap2_reduced)
  colnames(ohet)[1] <- "Heterozygosity" 
  colnames(ohet)[2] <- "Cluster"
  #make a boxplot
  box_plot <- boxplot(Heterozygosity ~ Cluster, data = ohet, notch = F, varwidth = TRUE)
}

#run the function
het_function(HAP1hwe, HAP2hwe, HEThwe)



