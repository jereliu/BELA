glmr_family <- 
  function(family_name = c("gaussian", "poisson", "binomial")){
    dist_family <- 
      paste0(family_name, "()") %>% 
      parse(text = .) %>% eval
    
    # log partition function and derivatives
    dist_family$partition <- partition(family_name)
    
    # TODO: add sufficient stat
    
    # return
    dist_family
  }

partition <- 
  function(family_name){
    if (family_name == "gaussian"){
      f = function(theta) theta^2/2
      d = function(theta) theta
      d2 = function(theta) 1 
    } else if (family_name == "poisson"){
      f = function(theta) theta
      d = function(theta) 1
      d2 = function(theta) 0 
    } else if (family_name == "binomial"){
      f = function(theta)  -log(1 - theta)
      d = function(theta)  1/(1 - theta) 
      d2 = function(theta)  1/(1 - theta)^2
    }
    
    out_list <- list(f = f, d = d, d2 = d2)
    
    return(out_list)
  }