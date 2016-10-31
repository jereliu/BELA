glmr_family <- 
  function(family_name = c("gaussian", "poisson", "binomial")){
    dist_family <- 
      paste0(family_name, "()") %>% 
      parse(text = .) %>% eval
    
    # log partition function and derivatives
    dist_family$partition <- partition(family_name)
    
    # TODO: add sufficient stat
    dist_family$sufficient <- sufficient(family_name)
    
    # negative log likelihood
    dist_family$negloglik <- 
      function(T, Theta, suffunc = dist_family$sufficient){
        A <- suffunc(Theta)
        sum(-Theta * T + A)
      }
    
    # return
    dist_family
  }

partition <- 
  function(family_name){
    if (family_name == "gaussian"){
      f = function(theta) theta^2/2
      d = function(theta) theta
      d2 = function(theta) theta*0 + 1 
    } else if (family_name == "poisson"){
      f = function(theta) exp(theta)
      d = function(theta) exp(theta)
      d2 = function(theta) exp(theta)
    } else if (family_name == "binomial"){
      f = function(theta)  log(1 + exp(theta))
      d = function(theta)  exp(theta)/(1 + exp(theta))
      d2 = function(theta)  exp(theta)/((1 + exp(theta))^2)
    }
    
    out_list <- list(f = f, d = d, d2 = d2)
    
    return(out_list)
  }

sufficient <-     # TODO: add sufficient stat
  function(family_name){
    suffunc <- function(Y) Y
    suffunc
  }