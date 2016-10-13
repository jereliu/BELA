using DataFrames
using ProgressMeter
cd("/Users/Jeremiah/GitHub/BELA/func/method/FISTA/julia")

####  1. Initialize Input Variables ####
# 1.0 setup task config
const link_type = "pos2"
const K = 100
const blow_factor = 1.1
const max_iter = 1000
const tol = 1e-5

# 1.1 read-in N
N = readtable("N.csv")
N = convert(Array, N)

n = size(N, 1)
p = size(N, 2)

# 1.2 initialize T and sigma
T = sum(N, 2)/10
sigma = ones(Int64, (size(N, 2), 1))
T_sigma = T * sigma'

# 1.3 initialize X and Y
# initialization using LU decomposition
X_init, Y_init = qr(sqrt(N./T_sigma))
Y_init = Y_init'

# 1.4 initialize Gamma_X, Gamma_Y
# identity case
Sigma_X = eye(n)
Sigma_Y = eye(p)

# 1.5 lambda
lambda = 1

#### 2. FISTA ####


#### 3. Experiments ####
Y = Y_init
X = X_init

D = compute_D(N, X, Y, T_sigma, link_type)
dL_X = L_deriv(X, Y, D, Sigma_inv_X, lambda, "X")
dL_Y = L_deriv(X, Y, D, Sigma_inv_Y, lambda, "Y")

Lip = init_Lip(1, "X", 1.,
               N, X, Y, T_sigma, Sigma_inv_X,
               dL_X, lambda, link_type)

blow_factor = 1.1

while true
  F_L_try =
    F_L(1, "X", Lip, N, X, Y, T_sigma,
    dL_X, lambda, link_type)
  Q_L_try =
    Q_L(1, "X", Lip, N, X, Y, T_sigma,
    dL_X, lambda, link_type)
  if F_L_try <= Q_L_try
    print("Lip found!, Lip = ", Lip)
    break
  end
  Lip = Lip * blow_factor
end
