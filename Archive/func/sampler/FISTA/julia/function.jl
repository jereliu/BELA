#### nabla_L, return derivative for entire X or Y
function compute_D(N, X, Y, T_sigma, link_type)
  # compute (non-factor specific) chain rule components for derivative of loss function
  # D = (T_sigma - N/f_Q) * f_Q_deriv
  Q = X*Y'
  # compute link function
  f_Q = link(Q, link_type)
  f_Q[f_Q .== 0] = 42.4242 # handle div by 0 case
  # compute link funciton derivative
  f_Q_deriv = link_deriv(Q, link_type)
  # compute D, a NxP matrix
  (T_sigma - N./f_Q) .* f_Q_deriv
end

function compute_D2(N, X, Y, T_sigma, link_type)
  # compute (non-factor specific) chain rule components for derivative of loss function
  # D = (T_sigma - N/f_Q) * f_Q_deriv
  Q = X*Y'
  # compute link function
  f_Q = link(Q, link_type)
  f_Q[f_Q .== 0] = 42.4242 # handle div by 0 case
  # compute link funciton derivative
  f_Q_deriv = link_deriv(Q, link_type)
  f_Q_deriv2 = link_deriv2(Q, link_type)

  # compute D, a NxP matrix
  N./(f_Q.^2) .* f_Q_deriv + (T_sigma - N./f_Q) .* f_Q_deriv2
end

function update_D()
  # TODO
end

function L_deriv(X, Y, D, Sigma_inv, lambda, wrt = "X")
  if wrt == "X"
    # dimension N x K
    D*Y + lambda * Sigma_inv * X
  elseif wrt == "Y"
    # dimension P x K
    D'*X + lambda * Sigma_inv * Y
  end
end

#### Q_L and F:
function F(fac_type,
          N, X, Y, T_sigma, Gamma,
          lambda, link_type)
  # TODO: not efficient. consider row-wise update
  Q = X*Y'
  # compute link function
  f_Q = link(Q, link_type)
  log_f_Q = log(f_Q)
  log_f_Q[log_f_Q.==-Inf] = 42.4242

  # compute loss
  if fac_type == "X"
    fac = X
  elseif fac_type == "Y"
    fac = Y
  else
    error("illegal fac_type")
  end

  sum(T_sigma .* f_Q - N .* log_f_Q) +
  lambda * vecnorm(Gamma * fac)^2
end

function F_full(N, X, Y, T_sigma, Gamma_X, Gamma_Y,
                lambda, link_type)
  # alternative version used to evaluate full objective
  Q = X*Y'
  # compute link function
  f_Q = link(Q, link_type)
  log_f_Q = log(f_Q)
  log_f_Q[log_f_Q.==-Inf] = 42.4242

  # compute loss
  sum(T_sigma .* f_Q - N .* log_f_Q) +
  lambda * (vecnorm(Gamma_X * X)^2 + vecnorm(Gamma_Y * Y)^2)
end

function F_L( # compute F(p_L) given L = Lip
  fac_idx, fac_type, Lip,
  N, X, Y, T_sigma, Gamma, dL,
  lambda, link_type)
  # fac_type: is the new factor X or Y
  # compute new factor
  if fac_type == "X"
    fac_old = X[fac_idx, :]
  elseif fac_type == "Y"
    fac_old = Y[fac_idx, :]
  else
    error("illegal fac_type")
  end

  fac_new =
    p_L(fac_old, fac_idx, Lip, dL)

  # update factor matrix then compute F_L
  X_new = copy(X)
  Y_new = copy(Y)
  if fac_type == "X"
    X_new[fac_idx, :] = fac_new
  elseif fac_type == "Y"
    Y_new[fac_idx, :] = fac_new
  end

  F(fac_type,
    N, X_new, Y_new, T_sigma, Gamma,
    lambda, link_type)
end

function Q_L(
  fac_idx, fac_type, Lip,
  N, X, Y, T_sigma, Gamma,
  dL, lambda, link_type)
  # fac_type: is the new factor X or Y
  # compute new and old factor
  if fac_type == "X"
    fac_old = X[fac_idx, :]
  elseif fac_type == "Y"
    fac_old = Y[fac_idx, :]
  else
    error("illegal fac_type")
  end

  fac_new =
    p_L(fac_old, fac_idx, Lip, dL)

  # compute Q_L
  fac_diff = fac_new - fac_old

  L_old = F(fac_type, N, X, Y, T_sigma, Gamma,
            lambda, link_type) # 1 x 1
  term_grad = (fac_diff * dL[fac_idx, :]')[1] # 1 x K * K x 1
  term_lips = 0.5 * Lip * norm(fac_diff)^2

  L_old + term_grad + term_lips
end


#### p_L: proximal operator given Lipschitz constant
function p_L(factor, i, Lip, dL)
  factor - dL[i, :]/Lip
end
