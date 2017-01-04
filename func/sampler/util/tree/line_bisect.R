solveLevelSet_bisect <- 
  function(
    # declare variable to solve
    var = "U", idx, 
    # current values
    U_cur, V_cur, T_suff, lambda,
    # directions
    U_dir, V_dir, 
    # functions and solutions
    slice_lik, target_lik)
  {
    # define function
    if (var == "U"){
      f <- function(a){
        U_add <- matrix(0, nrow = nrow(U_cur), ncol = ncol(U_cur))
        U_add[idx, ] <- U_dir
        target_lik(U_cur + a * U_add,
                   V_cur, T_suff, lambda) - slice_lik
      }
    } else {
      f <- function(a){
        V_add <- matrix(0, nrow = nrow(V_cur), ncol = ncol(V_cur))
        V_add[idx, ] <- V_dir
        target_lik(U_cur, V_cur + a * V_add, 
                   T_suff, lambda) - slice_lik
      }
    }
    
    scale_1 <- bisect(f, 0, 5)$root
    scale_2 <- bisect(f, 0, -5)$root
    
    # return two solutions
    if (var == "U"){
      solMat <- 
        rbind(U_cur[idx, ] + scale_1 * U_dir, 
              U_cur[idx, ] + scale_2 * U_dir)
    } else {
      solMat <- 
        rbind(V_cur[idx, ] + scale_1 * V_dir, 
              V_cur[idx, ] + scale_2 * V_dir)
    }
    solMat
  }

sampleLevelSet <- 
  function(
    # declare variable to solve
    var = "U", idx, 
    # current values
    U_cur, V_cur, T_suff, lambda,
    # directions
    U_dir, V_dir, 
    # functions and solutions
    slice_lik, target_lik, sample_scale = 5)
  {
    # define function
    if (var == "U"){
      f <- function(a){
        U_add <- matrix(0, nrow = nrow(U_cur), ncol = ncol(U_cur))
        U_add[idx, ] <- U_dir
        target_lik(U_cur + a * U_add,
                   V_cur, T_suff, lambda) - slice_lik
      }
    } else {
      f <- function(a){
        V_add <- matrix(0, nrow = nrow(V_cur), ncol = ncol(V_cur))
        V_add[idx, ] <- V_dir
        target_lik(U_cur, V_cur + a * V_add, 
                   T_suff, lambda) - slice_lik
      }
    }
    
    scl_l <- -sample_scale
    scl_u <- sample_scale
    a_prop <- scl_u
    
    while(f(a_prop) > 0){
      a_prop <- runif(1, scl_l, scl_u)
      if (a_prop > 0){
        scl_u <- a_prop
      } else {
        scl_l <- a_prop
      }
    }
    
    # return two solutions
    if (var == "U"){
      out <- U_cur[idx, ] + a_prop * U_dir
    } else {
      out <- V_cur[idx, ] + a_prop * V_dir
    }
    out
  }

sampleLevelSet_entire <- 
  function(
    # declare variable to solve
    var = "U",
    # current values
    U_cur, V_cur, T_suff, lambda,
    # directions
    U_dir, V_dir, 
    # functions and solutions
    slice_lik, target_lik, sample_scale = 5)
  {
    # define function
    if (var == "U"){
      f <- function(a){
        target_lik(U_cur + a * U_dir,
                   V_cur, T_suff, lambda) - slice_lik
      }
    } else {
      f <- function(a){
        target_lik(U_cur, V_cur + a * V_dir, 
                   T_suff, lambda) - slice_lik
      }
    }
    scl_l <- -sample_scale
    scl_u <- sample_scale
    a_prop <- scl_u
    
    while(f(a_prop) > 0){
      # if reject, modify sampling interval
      if (a_prop > 0){
        scl_u <- a_prop
      } else {
        scl_l <- a_prop
      }
      a_prop <- runif(1, scl_l, scl_u)
    }
    
    # return
    if (var == "U"){
      out <- U_cur + a_prop * U_dir
    } else {
      out <- V_cur + a_prop * V_dir
    }
    out
  }