##function for running LD from putative structural variants identified by lostruct##
##Barbara Langille March 2023

##this function uses chromosomes identified by lostruct as having putative structural variants
##this function makes a linkage disequilibrium heatmap

##you will need a .ld file of your entire chromosome - can do this in plink with following code:
##plink --file FILENAME  --r2 square --maf 0.01 --out FILENAME_LD --chr 1 #change chr number to the chromosome you are using





#packages needed
library(data.table)
library(superheat)





#Set the working directory for this script
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())

#this is the function - dont change anything unless you dont like the colors or the values for color ranges
LD_heatmap <- function(LDfile) {
  #read in the file
  CHR_LD  <- as.matrix(fread(LDfile, data.table = F))
  diag(CHR_LD) <- NA
  CHR_LD_Heatmap <- CHR_LD
  CHR_LD_Heatmap[upper.tri(CHR_LD_Heatmap, diag = T)] <- 0
  heatpal <-c("gray95", "blue", "purple", "red2", "red3")
  #make a plot
  superheat(CHR_LD_Heatmap, heat.pal = heatpal, heat.pal.values = c(0, 0.25, 0.5, 0.75, 1)) # may want to adjust these if your values are smaller
  
}


#run the function with your file
LD_heatmap(LDfile)

#If your plot has been working in the past and is now throwing the invalid graphics 
#state error, try resetting the graphics device by calling dev.off()
dev.off()


