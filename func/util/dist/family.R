glmr_family <- 
  function(family_name = c("gaussian", "poisson", "poisson_softplus", 
                           "poisson_reluappr", "binomial")){
    family_name <- match.arg(family_name)
    family_spec <- list()
    
    # check if want canonical parameterization instead of natural parameterization
    is_modify <- (length(grep("_", family_name)) > 0)
    if (is_modify){
      dist_family_name <- strsplit(family_name, "_")[[1]][1]
      link_name <- strsplit(family_name, "_")[[1]][2]
      family_spec <- c(family_spec, "link = 'identity'")
    } else {
      dist_family_name <- family_name
    }
    
    family_spec <- 
      do.call(paste, c(family_spec, sep = ",")) %>% paste0("(", ., ")")
    
    dist_family <- 
      paste0(dist_family_name, family_spec) %>% 
      parse(text = .) %>% eval
    
    if (is_modify){
      dist_family$linkinv <- 
        eval(parse(text = link_name))
    } else {
      # log partition function and derivatives
      dist_family$partition <- partition(dist_family_name)
      
      # TODO: add sufficient stat
      dist_family$sufficient <- sufficient(dist_family_name)
      
      # negative log likelihood
      dist_family$negloglik <- 
        function(T, Theta, partfunc = dist_family$partition$f){
          A <- partfunc(Theta)
          sum(-Theta * T + A)
        }
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