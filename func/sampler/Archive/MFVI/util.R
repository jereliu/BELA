conf_hypergeo_U <- 
  function (a, b, z, ip = 0) 
  {
    if (!is.complex(z)) 
      z = complex(real = z, imaginary = 0 * z)
    if (!is.complex(a)) 
      a = complex(real = a, imaginary = 0)
    if (!is.complex(b)) 
      b = complex(real = b, imaginary = 0)
    ans = 
      kummerM(z, a = a, b = b)/
      (cgamma(1 + a - b) * cgamma(1 - b)) + 
      (z^(1 - b)) * 
      kummerM(z, a = (1 + a - b), b = 2 - b)/
      (cgamma(a) * cgamma(b-1))
    ans
  }