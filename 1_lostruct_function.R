######################function for running lostruct#############################
##Barbara Langille March 2023

##this function assesses outlier loci using local PCAs 
##it then creates a file of the outliers on MDS1 and MDS2, and plots all windows for MDS axis 1

##this function uses a tped and map file generated in plink using the following call
##Do for each chromosome separately - change the --chr 1 call depending on which chromosome you want
##plink --file FILENAME --chr 1 --recode12 --transpose --out FILENAME_ch1 --allow-extra-chr

##this function uses windows=20 but you may want to change to 50 depending on how many snps you have
##I was using a 70K snp array




library(lfa) # this is for reading in the tped file
library(qpcR) # this is for concatonating windows with qpcR
library(lostruct) # this is the main program
library(ggplot2) # this is for making the plot
library(dplyr) # this is for the row_number function
library(data.table) # important for depreciating functions




#Set the working directory for this script
setwd("~/DIRECTORY")

#clear R's brain
rm(list=ls())

#this function is necessary for compiling all the outliers together - dont change anything
#Outliers function - get outliers in the dataset 
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

##this is the function - don't change anything in here!! 
#exception is if you want to change to a different window size 
lostruct_outliers <- function(tpedfile, MAPfile) {
  tped <- read.tped.recode(tpedfile, buffer.size = 5e+08) #read in the tped files
  #Compute the local PCA coordinates with 20 windows and 2 groups for each chromosome
  pcs <- eigen_windows(tped, win=20, k=2) #potentially want to change to win=50 if your snp file is bigger
  #Compute the distance matrix between windows
  pcdist <- pc_dist(pcs, npc=2)
  #MDS
  fit <- cmdscale(pcdist,eig=TRUE,k=2)
  #Get the map files
  MAP <- read.delim(MAPfile, stringsAsFactors = F, header = F, sep = "")
  #grab the start of all windows - in this case 20 - do for each chromosome
  MAP_Win20_start = MAP[seq(1, nrow(MAP), 20), ]
  #grab the end of all windows - in this case 20 - do for each chromosome
  MAP_Win20_end = MAP[seq(20, nrow(MAP), 20), ]
  #concatenate together
  MAP <- qpcR:::cbind.na(MAP_Win20_start, MAP_Win20_end, fit$points)
  #rename all columns
  colnames(MAP)[1] <- "CHRstart" 
  colnames(MAP)[2] <- "SNPstart"
  colnames(MAP)[3] <- "Xstart"
  colnames(MAP)[4] <- "POSstart"
  colnames(MAP)[5] <- "CHRend"
  colnames(MAP)[6] <- "SNPend"
  colnames(MAP)[7] <- "Xend"
  colnames(MAP)[8] <- "POSend"
  colnames(MAP)[9] <- "MDS1"
  colnames(MAP)[10] <- "MDS2"
  #IF your dataframe start and end are uneven, remove the last row
  MAP_nona <- MAP %>% filter(row_number() <= n()-1)
  #Write out the file with outliers and locations
  write.csv(MAP_nona, "lostruct_windows.csv")
  #Visualize MDS1
  n.sd <- 2 #outlier window cutoff is 2 or 3 standard deviations
  #get the cuttoff value -2 standard deviations away
  cutoff <- mean(MAP_nona$MDS1)+c(-1,1)*n.sd*sd(MAP_nona$MDS1)
  ##MDS plot with outliers against the window position
  ##red line in the plot is 2 standard deviations above and below the mean
  MDS1plot <- ggplot() + geom_point(data = MAP_nona, aes(x = POSstart, y = MDS1)) +
    geom_hline(yintercept = cutoff, linetype="dashed", color="red", size=0.5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Window position")+
    ylab(paste0("MDS1"))
  #saves as pdf but you can change to png if you prefer
  ggsave("MDS1plot.pdf", MDS1plot, width = 8, height = 6)
  #Outliers based on standard deviation from mean
  outlier_MSD1 <- outliers(MAP_nona$MDS1,2) #Using 2 standard deviations from mean
  outlier_MSD2 <- outliers(MAP_nona$MDS2,2) #Using 2 standard deviations from mean
  #Save outlier windows for MDS 1
  write.table(MAP_nona[which(abs(MAP_nona$MDS1) >= min(abs(outlier_MSD1))), ],
              file = "OutlierWindow_MDS1.txt",
              quote = F, row.names = F, col.names = T, sep="\t")
  #Save outlier windows for MDS 2
  write.table(MAP_nona[which(abs(MAP_nona$MDS2) >= min(abs(outlier_MSD2))), ],
              file = "OutlierWindow_MDS2.txt",
              quote = F, row.names = F, col.names = T, sep="\t")
  
}

# here is where you change the tped and map files that you are using
# will need to include path if you are not in the correct directory 
lostruct_outliers(tpedfile, MAPfile)




