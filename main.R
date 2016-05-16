require(magrittr)
require(dplyr)

source("./func/util/source_Dir.R")
sourceDir("./func")


#### 1. Generate Sample ####
N <- biom_sample(n_sample = 50, n_OTU = 20)

Q <- biom_sample(
  n_sample = 50, n_OTU = 20, 
  return = "Q")

#### 2. VI ####
df <- 
  apply(Q, 2, 
        function(x){
          x_temp <- x - min(x)
          x_temp/max(x_temp)
        }
  ) %>% t
out_str <- 
  RVineStructureSelect(df, familyset = 1)
out_cop <- 
  RVineTreePlot(data = df, RVM = out_str, method = "mle")
