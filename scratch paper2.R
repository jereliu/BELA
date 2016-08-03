#### 1 beta penalty ####
a_cat <- seq(0.01, 0.5, 0.01)
x <- seq(0, 1, 1e-3)
y_func <- function(alpha){
  beta <- 0.5 - alpha
  x <- seq(0, 1, 1e-3)
  y <- dbeta(x, alpha, beta)
  -y
}

#### 1.1 static ####

plot(x, y_func(0.01), type = "n",
     xlab = expression(sigma[i]), 
     ylab = expression(r(sigma[i]))
)

for (a in seq(0.01, 0.5, 0.01)){
  lines(x, y_func(a), lwd = 2, 
        col = rgb(0, 0, 0, 0.5))
}

plot(x, y_func(0.01), type = "n",
     xlab = expression(sigma[i]), 
     ylab = expression(r(sigma[i]))
)

for (a in seq(0.1, 0.2, 0.01)){
  lines(x, y_func(a), lwd = 2, 
        col = rgb(a*2, 0, 0))
}

#### 1.2 gif ####

for (a in a_cat){
  png(paste0("./temp_plot/beta_prior/beta_", a, ".png"))
  plot(x, y_func(0.01), type = "n",
       xlab = expression(sigma[i]), 
       ylab = expression(r(sigma[i])), 
       main = paste("sigma = ", a)
  )
  lines(x, y_func(a), lwd = 2)
  dev.off()
}


#### 2. objective function ####
n_cat <- c(0, 5, 20, 50, 100)
x_func <- function(n) seq(0, 200, 0.1)
y_func <- function(n){  
  x <- seq(0, 200, 0.1)
  y <- x - n * log(x) 
  y
}


plot(x_func(0), y_func(0), 
     xlim = c(0, 200), ylim = c(-500, 200),
     xlab = expression(z[ij]),
     ylab = expression(L(z[ij], D)),
     type = "n")

i <- 0
for (n in n_cat){
  i <- i + 1
  lines(x_func(n), y_func(n), col = i)
}
legend("bottomleft", lty = 1,
       col = 1:5, legend = paste0("n=", n_cat))

