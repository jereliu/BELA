find_V <- 
  function(U, Theta){
    # for Theta = UV^T, V = Theta^T * U * inv(U^TU)
    t(Theta) %*% U %*% inv(t(U) %*% U)
  }