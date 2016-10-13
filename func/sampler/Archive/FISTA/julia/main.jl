#### 1. main function ####
function FISTA(N, T, sigma, Sigma_X, Sigma_Y,
               Y_init, X_init,
               lambda = 1,
               link_type = "pos2",
               tol = 1e-5, blow_factor = 1.1,
               max_iter = 1000)
  #### 1. Initiatiations ####
  #### 1.1 Statistics ====
  # dimension information
  n = size(N, 1)
  p = size(N, 2)

  # condition factor
  T_sigma = T * sigma'
  obj_ideal = sum(N - N.*log(N./T_sigma + 1e-12))

  # prior covariance
  Sigma_inv_X = inv(Sigma_X)
  Sigma_inv_Y = inv(Sigma_Y)

  Gamma_X, _, _ = lu(Sigma_inv_X)
  Gamma_Y, _, _ = lu(Sigma_inv_Y)

  #### 1.2 Operation Container ====
  # Lipschitz Constant container
  Lip_X = ones(n)
  Lip_Y = ones(p)

  # Iteration Container
  X_cur = copy(X_init)
  X0_cur = copy(X_init)

  Y_cur = copy(Y_init)
  Y0_cur = copy(Y_init)

  obj_hist = zeros(max_iter)

  #### 2. Update version 1: FISTA ####
  obj_diff = Inf
  obj_cur = Inf
  fista_iter = 0
  t_cur = 1

  while abs(obj_diff) > tol
    # value iteration update
    fista_iter += 1
    obj_prev = copy(obj_cur)

    X0_prev = copy(X0_cur)
    Y0_prev = copy(Y0_cur)
    X_prev = copy(X_cur)
    Y_prev = copy(Y_cur)

    t_prev = t_cur

    # step size update
    t_cur = 0.5 * (1 + sqrt(1 + 4 * t_prev^2))

    #### 2.1 X ====
    # Update the Gradient Matrix
    D = compute_D(N, X_cur, Y_cur, T_sigma, link_type)
    dL_X = L_deriv(X_cur, Y_cur, D, Sigma_inv_X, lambda, "X")

    # Update X row-wise (i.e. by observation)
    #prog = Progress(n, 1)
    for i = 1:n
      # Update Lipschitz with Backtracking
      Lip_X[i] =
        find_Lip(i, "X", Lip_X[i], blow_factor,
                 N, X_cur, Y_cur, T_sigma,
                 Gamma_X, Sigma_inv_X, dL_X,
                 lambda, link_type)
      # Update x_i
      X0_cur[i,:] = p_L(X_cur[i,:], i, Lip_X[i], dL_X)

      # optionally, check objective
      #obj = F("X",
      #        N, X_cur, Y_cur, T_sigma, Gamma_X,
      #        lambda, link_type)
      #print(obj)
      #next!(prog)
    end

    X_cur = X0_cur +
           (X0_cur - X0_prev) * (t_prev - 1)/t_cur

    #### 2.2 Y ====
    # Update D, the Gradient Matrix component
    D = compute_D(N, X_cur, Y_cur, T_sigma, link_type)
    dL_Y = L_deriv(X_cur, Y_cur, D, Sigma_inv_Y, lambda, "Y")

    # Update Y row-wise (i.e. by OTU)
    #prog = Progress(p, 1)
    for j = 1:p
    # Update Lipschitz with Backtracking
      Lip_Y[j] =
        find_Lip(j, "Y", Lip_Y[j], blow_factor,
                 N, X_cur, Y_cur, T_sigma,
                 Gamma_Y, Sigma_inv_Y, dL_Y,
                 lambda, link_type)
      # Update y_j
      Y0_cur[j,:] = p_L(Y_cur[j,:], j, Lip_Y[j], dL_Y)

      # optionally, check objective
      #obj = F("X",
      #        N, X_cur, Y_cur, T_sigma, Gamma_X,
      #        lambda, link_type)
      #print(obj)
      #next!(prog)
    end

    Y_cur = Y0_cur +
           (Y0_cur - Y0_prev) * (t_prev - 1)/t_cur

    #### 2.3 Summarize ====
    obj_cur =
      F_full(N, X_cur, Y_cur, T_sigma,
             Gamma_X, Gamma_Y, lambda, link_type)

    obj_diff = obj_prev - obj_cur
    resid = N - link(X_cur * Y_cur', link_type).*T_sigma
    MSE = mean(sum(resid.^2))

    print("\nIter", fista_iter, ": ", obj_cur-obj_ideal,
          ". Error = ", MSE,
          ". Diff = ", round(obj_diff, 5))
    print("\n<-------------------------------->\n\n")

    if fista_iter > max_iter
      break
    end

    obj_hist[fista_iter] = obj_cur

    #### 2.4 Backstage dark magic ====
    #if obj_diff > 0
    #  # if direction is good, blow up the step size a bit
    #  Lip_X = Lip_X./blow_factor
    #  Lip_Y = Lip_Y./blow_factor
    #end
  end

  #### 3. Update version 2: Closed Form ####
  obj_diff = Inf
  obj_cur = Inf

  while abs(obj_diff) > tol
    # value iteration update
    obj_prev = copy(obj_cur)
    X_prev = copy(X_cur)
    Y_prev = copy(Y_cur)
    #### 2.1 X ====
    # Update the Gradient Matrix
    D = compute_D(N, X_cur, Y_cur, T_sigma, link_type)
    X_cur = -Sigma_X * D * Y_cur/lambda

    #### 2.2 Y ====
    # Update D, the Gradient Matrix component
    D = compute_D(N, X_cur, Y_cur, T_sigma, link_type)
    dL_Y = L_deriv(X_cur, Y_cur, D, Sigma_inv_Y, lambda, "Y")

    # Update Y row-wise (i.e. by OTU)
    #prog = Progress(p, 1)
    for j = 1:p
    # Update Lipschitz with Backtracking
      Lip_Y[j] =
        find_Lip(j, "Y", Lip_Y[j], blow_factor,
                 N, X_cur, Y_cur, T_sigma,
                 Gamma_Y, Sigma_inv_Y, dL_Y,
                 lambda, link_type)
      # Update y_j
      Y0_cur[j,:] = p_L(Y_cur[j,:], j, Lip_Y[j], dL_Y)

      # optionally, check objective
      #obj = F("X",
      #        N, X_cur, Y_cur, T_sigma, Gamma_X,
      #        lambda, link_type)
      #print(obj)
      #next!(prog)
    end

    Y_cur = Y0_cur +
           (Y0_cur - Y0_prev) * (t_prev - 1)/t_cur

    #### 2.3 Summarize ====
    obj_cur =
      F_full(N, X_cur, Y_cur, T_sigma,
             Gamma_X, Gamma_Y, lambda, link_type)

    obj_diff = obj_prev - obj_cur
    resid = N - link(X_cur * Y_cur', link_type).*T_sigma
    MSE = mean(sum(resid.^2))

    print("\nIter", fista_iter, ": ", obj_cur-obj_ideal,
          ". Error = ", MSE,
          ". Diff = ", round(obj_diff, 5))
    print("\n<-------------------------------->\n\n")

    if fista_iter > max_iter
      break
    end

    obj_hist[fista_iter] = obj_cur

    #### 2.4 Backstage dark magic ====
    #if obj_diff > 0
    #  # if direction is good, blow up the step size a bit
    #  Lip_X = Lip_X./blow_factor
    #  Lip_Y = Lip_Y./blow_factor
    #end
  end


  #### 4. Return ####

end


#### 2. Support function: Lipschitz Backtracking ####
# backtracking
function init_Lip(fac_idx, fac_type, Lip_init,
                  N, X, Y, T_sigma, Sigma_inv, dL,
                  lambda, link_type)
  ## guess1
  if fac_type == "X"
    fac_old = X[fac_idx, :]
  elseif fac_type == "Y"
    fac_old = Y[fac_idx, :]
  else
    error("illegal fac_type")
  end

  fac_new =
    p_L(fac_old, fac_idx, Lip_init, dL)

  # update factor matrix then compute F_L
  X_new = copy(X)
  Y_new = copy(Y)
  if fac_type == "X"
    X_new[fac_idx, :] = fac_new
  elseif fac_type == "Y"
    Y_new[fac_idx, :] = fac_new
  end

  # guess initial for Lip
  D_new = compute_D(N, X_new, Y_new, T_sigma, link_type)
  dL_new = L_deriv(X_new, Y_new, D_new, Sigma_inv, lambda, fac_type)

  guess1 = norm(dL_new - dL) / norm(fac_new - fac_old)

  ## guess2
  D2_new = compute_D2(N, X, Y, T_sigma, link_type)

  if fac_type == "X"
    guess2 = D2_new[fac_idx, :] * sum(Y.^2, 2)
  elseif fac_type == "Y"
    guess2 = D2_new[:, fac_idx]' * sum(X.^2, 2)
  else
    error("illegal fac_type")
  end

  guess2 = guess2[1]

  #print("\n\nguess1 = ", guess1, ", guess2 = ", guess2, "\n")
  ## max
  max(guess1, guess2, 1)
end


function find_Lip(fac_idx, fac_type, Lip_init, blow_factor,
                  N, X, Y, T_sigma, Gamma, Sigma_inv,
                  dL, lambda, link_type)
  # blow_factor: factor to blow up for next Lipschits constant search
  # initialize
  if Lip_init == 1
    Lip = init_Lip(fac_idx, fac_type, Lip_init,
                   N, X, Y, T_sigma, Sigma_inv,
                   dL, lambda, link_type)
  else
    Lip = Lip_init
  end

  # backtracking for proper Lipschitz constant
  while true
    F_L_try =
      F_L(fac_idx, fac_type, Lip,
          N, X, Y, T_sigma, Gamma, dL,
          lambda, link_type)
    Q_L_try =
      Q_L(fac_idx, fac_type, Lip,
          N, X, Y, T_sigma, Gamma,
          dL, lambda, link_type)
    if F_L_try <= Q_L_try
      #print("Lip found!, Lip = ", Lip)
      break
    end
    Lip = Lip * blow_factor
  end

  Lip
end
