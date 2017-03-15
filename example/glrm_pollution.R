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
library(NMF)
library(Rtsne)
source("./func/util/source_Dir.R")
sourceDir("./func")


addr_targ <- 
  paste0("../../Dropbox (Personal)/Research/Harvard/Lorenzo/1. MComp/Paper/plot/")

#### 1. Data Import ####
dStudy <- 
  fread("./func/data/real/airpollution.csv", 
        showProgress = FALSE, data.table = FALSE)

# extract names
elem_name <- 
  c("Ca","Cu", "Cl", "Fe", "Pb", "Br",
    "K", "Na", "Mg", "Mn", "Ni", "V", 
    "Si", "Zn", "Ti", "no2")
namList <- paste0("log", elem_name, "Ratio_out")

Y <- na.omit(dStudy[grep("MAD-001", dStudy$SiteID), namList])

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
    do.call(data.frame,lapply(Y, function(x) replace(x, is.infinite(x), NA))) %>% 
    na.omit()
  rm_idx <- which(Y_plot[, 1] > 0)
  if (length(rm_idx) > 0) Y_plot <- Y_plot[-rm_idx, ]
  
  png(paste0(addr_targ, "corr.png"), 960, 960)
  pairs.panels(Y_plot)
  dev.off()
  
  my_palette <- 
    colorRampPalette(c("white", "yellow", "red"))(n = 299)
  
  col_breaks <- 
    c(seq(0, 0.01,length=100),  # for red
      seq(0.011,0.1,length=100),           # for yellow
      seq(0.101,1,length=100))             # for green
  
  pdf(paste0(addr_targ, "heatmaps.pdf"), 
      width = 8, height = 8)
  heatmap.2(cor(Y_plot),
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

#### 1. Baseline Model ####
V_est <- factanal(Y, 7)$loadings
class(V_est) <- "matrix"
V_est[V_est < 0] <- 0

#### 2. Naive Model ####
Y_plot <- 
  do.call(data.frame,lapply(Y, function(x) replace(x, is.infinite(x),NA))) %>% 
  na.omit()
# rm_idx <- which(Y_plot[, 1] > 0)
Y <- na.omit(Y_plot) %>% as.matrix()
Y_ana <- exp(Y) %>% scale(center = FALSE)

rec_pollution <-
  glrm(Y_ana, 
       lambda = 1, 
       k = 5, 
       true_par = NULL,
       init = NULL, 
       samplr_name = "hmc_stan",
       family_name = "gaussian",
       prior_name = "dirichlet",
       iter_max = c(1e3, 1e4),
       record_freq = 1
  )

error <- 
  apply(rec_pollution$Theta, 1, 
        function(Theta) 
          mean((Y_ana - log(Theta))^2/Y_ana^2)
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
idx <- 
  round((dim(rec_pollution$V)[1] * 0.9):(dim(rec_pollution$V)[1]))

V_est_vm <- 
  lapply(idx, function(id){
    out <- varimax(rec_pollution$V[id, ,])$loadings
    class(out) <- "matrix"
    out <- t(t(abs(out))/colSums(abs(out)))
    t(out * (out > 0.1))
  }) %>% do.call(rbind, .)

V_est_clust <- kmeans(V_est_vm, 8)
# clusplot(V_est_vm, 
#          V_est_clust$cluster, color=TRUE, shade=TRUE, 
#          labels=2, lines=0)

out <- aggregate(V_est_vm, by=list(V_est_clust$cluster),FUN = mean)
out <- t(out[1:8, -1])
rownames(out) <- elem_name



#### 1.2 for dirichlet ----
idx <- 
  round((dim(rec_pollution$V)[1] * 0.9):(dim(rec_pollution$V)[1]))

V_est_all <- 
  lapply(idx, function(i) lof(rec_pollution$V[i, , ])) %>% 
  simplify2array %>% apply(1, rbind) %>% unique

image(lof(t(V_est_all)))

V_est_clust <- Rtsne(V_est_all, dims = 2)
plot(V_est_clust$Y)

rownames(out) <- elem_name


#### 3. Plot ====

# before rotation 
csv <- V_est
# V_est_vm <- varimax(t(V_est))$loadings
# csv <- t(t(abs(V_est_vm))/colSums(abs(V_est_vm)))
class(csv) <- "matrix"
rownames(csv) <- elem_name
colnames(csv) <- 
  c("Traffic Exhaust/Road Dust", 
    "Metal Processing Plant", 
    "Road Salt", 
    "Biomass Burning", 
    "Secondary Road Exhaust", 
    "Fertilizer Factory",
    "Residual Oil Combustion")

csv.m <- melt(csv, id.vars="V1")
csv.m$Var2 <- as.character(csv.m$Var2)

qplot(x=Var1, y=Var2, data=csv.m, 
      fill=value, geom="tile") + 
  theme(axis.text.x = 
          element_text(angle = 90, hjust = 1)) + 
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red") + 
  labs(x = "Source", y = "Chemicals",
       title = "Latent Pollution Source, Before Rotation")

ggsave(paste0(addr_targ, "ex1_pmf.pdf"))

# after rotation ====
csv <- apply(out, 2, function(x) x/sum(x)) 
class(csv) <- "matrix"
colnames(csv) <- 
  c("Motor Vehicles", "Residual oil combustion", 
    "Metal Plant", "4", "5", 
    "6", "7", "8")
rownames(csv) <- elem_name

csv.m <- melt(t(csv), id.vars="V1")
csv.m$Var1 <- as.character(csv.m$Var1)

qplot(y=Var1, x=Var2, data=csv.m, 
      fill=value, geom="tile") + 
  theme(axis.text.x = 
          element_text(angle = 90, hjust = 1)) + 
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red") + 
  ylab(10:2) + 
  labs(y = "Pollution Source", x = "Chemical",
       title = "Source Emission Profile, After Rotation")

ggsave(paste0(addr_targ, "polfactor_afterot.pdf"))
