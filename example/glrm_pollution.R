rm(list = ls()) #WoG

library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
library(gplots)
library(data.table)
library(psych)

library(reshape2)
library(ggplot2)
library(cluster)
library(mclust)
source("./func/util/source_Dir.R")
sourceDir("./func")


addr_targ <- 
  paste0("../../Dropbox (Personal)/Research/Harvard/Lorenzo/",
         "1. BayesOpt/Report/Progress/", 
         "2017_Nov_Week_4/plot/")

#### 1. Data Import ####
dStudy <- 
  fread("./func/data/real/airpollution.csv", 
        showProgress = FALSE, data.table = FALSE)

# extract names
elem_name <- 
  c("S", "K", "Ca", "Fe", "Zn", "Cu", "Ti", "Al", "Pb", "Cl", "V", "Ni")
namList <- paste0("log", elem_name, "Ratio")

Y <- dStudy[, namList]

diag_plot <- FALSE
if (diag_plot){
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
}

#### 2. Naive Model ####
Y_plot <- 
  do.call(data.frame,lapply(Y, function(x) replace(x, is.infinite(x),NA)))
rm_idx <- which(Y_plot[, 1] > 0)
Y <- na.omit(Y_plot[-rm_idx, ]) %>% as.matrix()
Y_ana <- Y[sample(nrow(Y), 200), ]


rec_pollution <- 
  glrm(Y_ana, 
       lambda = 10, 
       k = 10, 
       true_par = NULL,
       init = NULL, 
       init_MAP = TRUE,
       samplr_name = "hmc_stan",
       family_name = "gaussian",
       iter_max = c(1e3, 1e4),
       record_freq = 1
  )
save(rec_pollution, 
     file = "./result/pollution/rec_pollution_naive.RData")

# prediction
load("./result/pollution/rec_pollution_naive.RData")
pred_num <- 100
max_iter <- dim(rec_pollution$U)[1]
pred_iter <- seq(max_iter/pred_num, max_iter, length.out = pred_num)

pred_Theta <- 
  lapply(pred_iter, 
         function(i){
           apply(rec_pollution$Theta[1:i, , ], 
                 c(2, 3), mean)
         }
  ) %>% abind(along = 3)

pred_error <- 
  apply(pred_Theta, 3, 
        function(theta)
          sqrt(sum((theta - rec_pollution$Y)^2)/
                 sum((rec_pollution$Y - mean(rec_pollution$Y))^2))
  )

pdf(paste0(addr_targ, "polfactor_prediction.pdf"), 
    height = 4, width = 10)
par(mfrow = c(1, 2))
plot(rec_pollution$lp__[20:length(rec_pollution$lp__)], 
     type = "l", 
     xlab = "Iteration", 
     ylab = "log posterior likelihood")

plot(pred_iter[-1], pred_error[-1], 
     xlab = "Iterations", ylab = "Standardized RMSE",
     type = "l", main = "Prediction Error")
dev.off()

# 
idx <- (dim(rec_pollution$V)[1] * 0.8):(dim(rec_pollution$V)[1])

V_est_vm <- 
  lapply(idx, function(id){
    out <- varimax(rec_pollution$V[id, ,])$loadings
    class(out) <- "matrix"
    out <- t(t(abs(out))/colSums(abs(out)))
    t(out * (out > 0.1))
  }) %>% do.call(rbind, .)

V_est_clust <- kmeans(V_est_vm, 40)
# clusplot(V_est_vm, 
#          V_est_clust$cluster, color=TRUE, shade=TRUE, 
#          labels=2, lines=0)

out <- aggregate(V_est_vm, by=list(V_est_clust$cluster),FUN = mean)
out <- t(out[1:10, -1])
rownames(out) <- rev(elem_name)



#### 3. Plot ====

# before rotation 
csv <- t(t(abs(V_est))/colSums(abs(V_est)))
rownames(csv) <- elem_name

csv.m <- melt(t(csv), id.vars="V1")
csv.m$Var2 <- as.character(csv.m$Var2)

qplot(x=Var1, y=Var2, data=csv.m, 
      fill=value, geom="tile") + 
  theme(axis.text.x = 
          element_text(angle = 90, hjust = 1)) + 
  labs(x = "Factors") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red") + 
  labs(x = "Factors", y = "Pollutant",
       title = "Latent Pollution Source, Before Rotation")

ggsave(paste0(addr_targ, "polfactor_beforot.pdf"))

# after rotation ====
csv <- t(t(out)/colSums(out))
class(csv) <- "matrix"
rownames(csv) <- rev(elem_name)

csv.m <- melt(t(csv), id.vars="V1")
csv.m$Var1 <- as.character(csv.m$Var1)

qplot(x=Var1, y=Var2, data=csv.m, 
      fill=value, geom="tile") + 
  theme(axis.text.x = 
          element_text(angle = 90, hjust = 1)) + 
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red") + 
  scale_x_discrete(labels=1:10) + 
  labs(x = "Factors", y = "Pollutant",
       title = "Latent Pollution Sources (V), After Rotation, Scaled")

ggsave(paste0(addr_targ, "polfactor_afterot.pdf"))
