# adapted from "GenerateCliquesCovariance" function 
cov_sparse <- 
  function (ncliques, cliquesize, theta) 
  {
    p <- ncliques * cliquesize
    sizes <- rep(cliquesize, ncliques)
    Sigma <- matrix(0, p, p)
    lims <- c(0, cumsum(sizes))
    for (i in seq(ncliques)) {
      ii <- (lims[i] + 1):lims[i + 1]
      signs <- 2 * (rnorm(sizes[i] * (sizes[i] - 1)/2) < 0) - 
        1
      Sigma[ii, ii][upper.tri(Sigma[ii, ii])] <- signs * theta
    }
    Sigma <- Sigma + t(Sigma)
    eig <- eigen(Sigma, symmetric = T)
    shift <- (max(eig$val) - p * min(eig$val))/(p - 1)
    #cat("Shifting eigenvalues by ", shift, fill = T)
    diag(Sigma) <- diag(Sigma) + shift
    A <- eig$vect %*% diag(sqrt(eig$val + shift)) %*% t(eig$vect)
    Sigma = Sigma
  }