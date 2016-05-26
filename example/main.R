require(magrittr)
require(dplyr)

source("./func/util/source_Dir.R")
sourceDir("./func")

#### 1. Simple Copula Update ####
N <- 
  biom_sample(n_sample = 50, n_OTU = 20)

Q <- 
  biom_sample(n_sample = 50, n_OTU = 20, 
              return = "Q")
