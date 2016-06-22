pred_dist <- 
  function(info, N){
    pred_mat <- 
      (info$mean$sigma * t(info$mean$Q2)) %>% 
      t %>% apply(1, function(x) x/sum(x))
    orig_mat <- 
      apply(N, 1, function(x) x/sum(x))
    orig_mat_norm <- 
      orig_mat %>% abs %>% mean
    
    (orig_mat-pred_mat) %>% 
      divide_by(orig_mat_norm) %>% abs %>% mean
  }


ELBO <- function(lambda, prior, info, 
                 par_name = "T",verbose = FALSE){
  
}