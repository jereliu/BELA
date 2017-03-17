# lof form
lof <- function(V, tol = 5e-3){
  dim2 <- ncol(V)
  src_score <- 
    apply(V, 2, 
          function(v) paste(as.numeric(v>tol), collapse = "")) 
  src_order <- src_score %>% order %>% rev
  V[, src_order]
}