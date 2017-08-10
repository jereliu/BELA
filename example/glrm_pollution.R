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
  paste0("../../Dropbox (Personal)/Research/",
         "Harvard/Thesis/Lorenzo/1. MComp/Paper/plot/")

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
  do.call(data.frame,
          lapply(Y, function(x) 
            replace(x, is.infinite(x),NA))) %>% 
  na.omit()
# rm_idx <- which(Y_plot[, 1] > 0)
Y <- na.omit(Y_plot) %>% as.matrix()
Y_ana <- exp(Y) %>% scale(center = FALSE)
n_iter <- 5e2

# dirichlet model
rec_pollution <-
  glrm(Y_ana, 
       lambda = 1, 
       k = 10, 
       true_par = NULL,
       init = NULL, 
       samplr_name = "hmc_stan",
       family_name = "gaussian",
       prior_name = "dirichlet_sparse",
       iter_max = c(1e3, n_iter),
       record_freq = 10
  )

rec_pollution_gam <-
  glrm(Y_ana, 
       lambda = 1, 
       k = 10, 
       true_par = NULL,
       init = NULL, 
       samplr_name = "hmc_stan",
       family_name = "gaussian",
       prior_name = "dirichlet_sparse_gam_sample",
       iter_max = c(1e3, n_iter),
       record_freq = 10
  )

# gaussian model
# rec_pollution <-
#   glrm(scale(Y),
#        lambda = 1, 
#        k = 5, 
#        true_par = NULL,
#        init = NULL, 
#        samplr_name = "hmc_stan",
#        family_name = "gaussian",
#        prior_name = "sparse",
#        iter_max = c(1e3, n_iter),
#        record_freq = 1
#   )

error <-
  apply(rec_pollution$Theta, 1, 
        function(Theta) 
          mean((Y_ana - log(Theta))^2/Y_ana^2)
  )


V_trace <- 
  apply(rec_pollution$V, 1, function(X) sum(X^2))
V_vol <- 
  apply(rec_pollution$V, 1, function(V) determinant(t(V) %*% V)$modulus)

rec_pollution$V[5e2, , ] %>% round(3)
rec_pollution$V[1e3, , ] %>% round(3)

V_trace <- 
  apply(rec_pollution$V[n_iter/2:n_iter, ,], 1, 
        function(V) sum(V[, 5]^2))

save(rec_pollution, 
     file = "./result/pollution/rec_pollution_naive.RData")

save(rec_pollution_gam, 
     file = "./result/pollution/rec_pollution_gam.RData")


# prediction
load("./result/pollution/rec_pollution_naive.RData")
load("./result/pollution/rec_pollution_gam.RData")

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

# varimax
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

# prediction
factor_mat <- list()
for (k in 1:10){
  factor_mat[[k]] <- list()
  for (iter in 1:n_iter){
    factor_mat[[k]][[iter]] <- 
      rec_pollution$U[iter+1, , k] %*% 
      t(rec_pollution$V[iter+1, , k])
  }
  factor_mat[[k]] <- 
    abind(factor_mat[[k]], along = 3)
}

scale2 <- 
  function(x) (x-min(x))/(max(x) - min(x))

plot(scale2(Y_ana[, 1]), 
     col = 2, type = "l", lwd = 2)

for (p in 1:16){
  lines(factor_mat[[3]][,p,5e2] %>% scale2,
        col = p)
}

#### 1.2 check population distribution of dirichlet factors ----
idx <-
  round((dim(rec_pollution_gam$V)[1] * 0.9):
          (dim(rec_pollution_gam$V)[1]))

# fixed effect prediction
time <- 1:nrow(Y)
X <- cbind(1, ns(time, df = 2))
B_est <- apply(rec_pollution_gam$B, c(2,3), mean)
Yh_t <- X %*% t(B_est)

for (i in 1:ncol(Y_ana)){
  plot(Y_ana[, i], 
       type = "l", col = 2, main = elem_name[i])
  lines(Yh_t)
}

# overall factor distribution
V_est_all <- 
  lapply(idx, function(i) lof(rec_pollution$V[i, , ])) %>% 
  simplify2array %>% apply(1, rbind) %>% unique

image(lof(t(V_est_all)), col = rev(heat.colors(12)))

V_est_clust <- Rtsne(V_est_all, dims = 2)
plot(V_est_clust$Y)

rownames(out) <- elem_name

# factor distribution by factor index
rec_pollution_gam <- rec_pollution 
K <- dim(rec_pollution_gam$U)[3]
idx <- 1:(dim(rec_pollution_gam$V)[1]/2 - 3)
# round((dim(rec_pollution_gam$V)[1] * 0.9):
#         (dim(rec_pollution_gam$V)[1]))

factor_imp <- 
  apply(rec_pollution_gam$gamma_rank_cumprod[idx, ,], 
        2, mean)

V_est_k <- 
  lapply(
    1:K, function(k) 
      t(rec_pollution_gam$V[idx, , k])
  )

for (k in 1:K){
  #pdf(paste0(addr_targ, "factor_", k, ".pdf"),
  #    width = 10, height = 5)
  #par(cex = 1.3)
  image(V_est_k[[k]], 
        col = rev(heat.colors(12)), 
        main = paste0("Factor ", k, 
                      ", Importance ", round(factor_imp[k], 4)),
        xlab = "Element Probability",
        ylab = "Iteration",
        axes = FALSE)
  axis(1, at = seq(0, 1, length.out = length(elem_name)), 
       labels= elem_name)
  axis(2, at = seq(0, 1, length.out = 5), 
       labels= seq(0, (dim(rec_pollution_gam$U)[1]-1)*10, 
                   length.out = 5))
  #dev.off()
}



#### 3. Plot ====

# before rotation 
csv <- V_est <- apply(rec_pollution$V, c(2,3) , mean)

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
