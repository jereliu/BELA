rm(list = ls()) #WoG

library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
library(gplots)
library(data.table)
source("./func/util/source_Dir.R")
sourceDir("./func")

addr_targ <- 
  paste0("../../Dropbox (Personal)/", 
         "Research/Harvard/Lorenzo/", 
         "1. BayesOpt/Report/7. Alg Review/plot/")

#### 1. Data Import ####
dStudy <- 
  fread("./func/data/real/airpollution.csv", 
        showProgress = FALSE, data.table = FALSE)

# extract names
elem_name <- 
  c("S", "K", "Ca", "Fe", "Zn", "Cu", "Ti", "Al", "Pb", "Cl", "V", "Ni")
namList <- paste0("log", elem_name, "Ratio")

Y <- dStudy[, namList]

# density
pdf(paste0(addr_targ, "pollution_histogram.pdf"), 
    width = 12, height = 8)
par(mfrow = c(3, 4))
for (elem in namList){
  #plot(density(log(Y[, elem])), main = paste0("log", elem))
  hist(Y[, elem], main = paste0("log", elem, "Ratio"))
}
dev.off()

# correlation
Y_plot <- 
  do.call(data.frame,lapply(Y, function(x) replace(x, is.infinite(x),NA)))
rm_idx <- which(Y_plot[, 1] > 0)

png(paste0(addr_targ, "corr.png"), 960, 960)
pairs.panels(Y_plot[-rm_idx, ])
dev.off()

my_palette <- 
  colorRampPalette(c("white", "yellow", "red"))(n = 299)

col_breaks <- 
  c(seq(0, 0.01,length=100),  # for red
    seq(0.011,0.1,length=100),           # for yellow
    seq(0.101,1,length=100))             # for green

pdf(paste0(addr_targ, "heatmaps.pdf"), 
    width = 8, height = 8)
heatmap.2(cor(Y),
          cellnote = round(cor(Y), 3),  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          Rowv = FALSE, Colv = FALSE,
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks    # enable color transition at specified limits
)
dev.off()

# censoring 
apply(dStudy[, elem_name], 2, function(x) sum(x==0)/length(x))*100
