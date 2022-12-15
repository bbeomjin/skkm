kernelMatrix = function(x, y, kernel = "gaussian", kparam = 1.0)
{
  kernel = match.arg(kernel, c("linear", "poly", "gaussian", "spline", "anova_gaussian") )
  x = as.matrix(x)
  y = as.matrix(y)
  p = ncol(x)
  
  if (ncol(x) == 0) {
    x = matrix(0, nrow = nrow(x), ncol = 1)
  }
  
  if (ncol(y) == 0) {
    y = matrix(0, nrow = nrow(y), ncol = 1)
  }
  
  if (kernel == "poly") {
    K = (x %*% t(y) + 1.0)^kparam
  } else if(kernel == "gaussian") {
    normx = rowSums(x^2)
    normy = rowSums(y^2)
    temp = x %*% t(y)
    temp = (-2.0 * temp) + outer(normx, rep(1.0, nrow(y)), "*") + outer(rep(1.0, nrow(x)), normy, "*")
    K = exp(-temp * kparam)
    # K = kernlab:::kernelMatrix(rbfdot(sigma = kparam), x, y)
  } else if (kernel == "linear") {
    # K = tcrossprod(x, y)
    K = x %*% t(y)
  } else if (kernel == "anova_gaussian") {
    K = 0
    for (d in 1:p) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, kernel = "gaussian", kparam = kparam)
      K = K + K_temp
    }
  } else if (kernel == "spline") {
    if (ncol(x) > 1 | ncol(y) > 1) {
      error("Spline kernel allows only one dimensional data")
    } else {
      K1x = (x - 1 / 2)
      K1y = (y - 1 / 2)
      K2x = (K1x^2 - 1 / 12) / 2
      K2y = (K1y^2 - 1 / 12) / 2
      ax = x %x% matrix(1, 1, nrow(y))
      ay = y %x% matrix(1, 1, nrow(x))
      b = abs(ax - t(ay))
      K1 = K1x %x% t(K1y)
      K2 = K2x %x% t(K2y) - ((b - 1 / 2)^4 - (b - 1 / 2)^2 / 2 + 7 / 240) / 24
      K = list(K1 = K1, K2 = K2)
    } 
  } else {
    K = NULL
  }
  return(K)
}



anovaKernel.gaussian = function(x, y, kernel, kparam)
{
  out = list()
  kernel = match.arg(kernel, c("gaussian", "gaussian-2way"))
  x = as.matrix(x)
  y = as.matrix(y)
  dimx = ncol(x)

  anovaKernel = lapply(1:dimx, function(j) {
                        kernelMatrix(x[, j, drop = FALSE], y[, j, drop = FALSE], 
                                     kernel = "gaussian", kparam = kparam)
                      })
  names(anovaKernel) = paste0("x", 1:dimx)
  numK = dimx
  
  if (grepl("2way", kernel)) {
    nint = dimx * (dimx - 1) / 2
    anovaKernel_int = vector(mode = "list", nint)
    index = 0
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anovaKernel[[i]]
        B = anovaKernel[[j]]
        anovaKernel_int[[index]] = A * B
      }
    }
    comb_name = combn(names(anovaKernel), 2)
    names(anovaKernel_int) = apply(comb_name, 2, paste0, collapse = "")
    anovaKernel = c(anovaKernel, anovaKernel_int)    
    numK = numK + nint
  }
  out$K = anovaKernel
  out$numK = numK
  return(out)
}

anovaKernel.linear = function(x, y, kernel, kparam = NULL)
{
  out = list()
  kernel = match.arg(kernel, c("linear"))
  x = as.matrix(x)
  y = as.matrix(y)
  dimx = ncol(x)

  anovaKernel = lapply(1:dimx, function(j) {
                        kernelMatrix(x[, j, drop = FALSE], y[, j, drop = FALSE], 
                                     kernel = "linear", kparam = NULL)
                      })
  names(anovaKernel) = paste0("x", 1:dimx)
  out$K = anovaKernel
  out$numK = dimx
  return(out)
}

anovaKernel.poly = function(x, y, kernel, kparam = 1)
{
  out = list()
  kernel = match.arg(kernel, c("poly"))
  x = as.matrix(x)
  y = as.matrix(y)
  dimx = ncol(x)

  anovaKernel = lapply(1:dimx, function(j) {
                        kernelMatrix(x[, j, drop = FALSE], y[, j, drop = FALSE], 
                                     kernel = "poly", kparam = kparam)
                      })
  names(anovaKernel) = paste0("x", 1:dimx)
  out$K = anovaKernel
  out$numK = dimx
  return(out)
}


anovaKernel.spline = function(x, y, kernel, kparam)
{
  out = list()
  kernel = match.arg(kernel, c("spline", "spline-t", "spline-2way", "spline-t-2way"))
  x = as.matrix(x)
  y = as.matrix(y)
  dimx = ncol(x)
  
  if (kernel == "spline") {
    numK = 2 * dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, kernel = "spline")
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep="")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep="")
    }
    
  } else if (kernel == 'spline-2way') {
    numK = (2 * dimx) + (2 * dimx * (2 * dimx - 1) / 2 - dimx)
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    # main effects
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, kernel = "spline")
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep = "")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep = "")
    }
    # two-way interactions
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A_linear = as.matrix(anova_kernel[[2 * i - 1]])
        A_smooth = as.matrix(anova_kernel[[2 * i]])
        B_linear = as.matrix(anova_kernel[[2 * j - 1]])
        B_smooth = as.matrix(anova_kernel[[2 * j]])
        anova_kernel[[index]] = A_linear * B_linear
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_linear * B_smooth
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " smooth", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_linear
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_smooth
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " smooth", sep = "")
      }
    }
  } else if (kernel == "spline-t") {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, kernel = "spline")
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
  } else if (kernel == 'spline-t-2way') {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, kernel = "spline")
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
    
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  }
  out$K = anova_kernel
  out$numK = numK 
  return(list(K = anova_kernel, numK = numK))
}


combine_kernel = function(anovaKernel, 
                          theta = rep(1 / sqrt(anovaKernel$numK),
                                      anovaKernel$numK))
{
  K = 0
  for (v in 1:length(theta)) {
    K = (K + theta[v] * anovaKernel$K[[v]])
  }
  return(K)
}

# make_a_vec = function(anovaKernel, clusters) {
#   uc = sort(unique(clusters))
#   d = anovaKernel$numK
#   n = length(clusters)
#   K1 = 2 * n * sapply(anovaKernel$K, function(x) sum(diag(x)))
#   K2 = 2 * sapply(anovaKernel$K, sum)
#   
#   K3 = K4 = numeric(d)
#   for (j in 1:d) {
#     r_part = lapply(uc, function(i) {
#       gind = clusters == i
#       temp_K = anovaKernel$K[[j]][gind, gind]
#       temp_K3 = sum(gind) * sum(diag(temp_K))
#       temp_K4 = sum(temp_K)
#       return(list(temp_K3, temp_K4))
#     })
#     K3[j] = sum(sapply(r_part, "[[", 1))
#     K4[j] = sum(sapply(r_part, "[[", 2))
#   }
#   return(K1 - K2 - K3 + K4)
# }

GetWCD = function(anovaKernel, clusters, weights = NULL)
{
  uc = unique(clusters)
  d = anovaKernel$numK
  wcd = numeric(d)
  for (v in 1:d) {
    K = anovaKernel$K[[v]]
    
    for (g in 1:length(uc)) {
      ind = clusters == uc[g]
      swt = weights[ind]
      subK = K[ind, ind]
      # wcd[v] = wcd[v] + sum(diag(subK)) - (1 / sum(swt)) * sum((subK * tcrossprod(swt)))
      wcd[v] = wcd[v] + sum(swt * diag(subK)) - (1 / sum(swt)) * sum((subK * tcrossprod(swt)))
    }
  }
  return(wcd)
  # return(list(td = td, wcd = wcd))
}


updateCs = function(K, clusters, weights, maxiter = 100) {
  
  # Initialization
  clusters0 = clusters
  
  for (i in 1:maxiter) {
    uc = unique(clusters0)
    RKHS_dist = sapply(uc, function(g) {
      
      ind = clusters0 == g
      swt = weights[ind]
      
      Kxx = diag(K)
      
      # Kxy = list()
      # Kxy$K = lapply(anovaKernel$K, function(x) x[ind, , drop = FALSE] * swt)
      # Kxy = colSums(combine_kernel(Kxy, theta))
      Kxy = colSums(K[ind, , drop = FALSE] * swt)
      # Kxy = (K[, ind, drop = FALSE] %*% swt) / sum(swt)
      
      # Kyy = list()
      # Kyy$K = lapply(anovaKernel$K, function(x) x[ind, ind, drop = FALSE] * tcrossprod(swt))
      # Kyy = sum(combine_kernel(Kyy, theta))
      Kyy = sum(K[ind, ind, drop = FALSE] * tcrossprod(swt))
      # Kyy = sum(drop(crossprod(K[ind, ind], swt)) * swt) / sum(swt)^2
      
      return(Kxx - (2 * Kxy / sum(swt)) + (Kyy / sum(swt)^2))
    })
    clusters = uc[apply(RKHS_dist, 1, which.min)]
    if (sum(clusters != clusters0) == 0) {
      break
    } else {
      clusters0 = clusters
    }
  }
  return(list(clusters = clusters, iteration = i, RKHS_dist = RKHS_dist))
}


soft_threshold = function(x, delta) {
  w = sign(x) * pmax(abs(x) - delta, 0)
  return(w)
}

normalization = function(x) {
  return(x / sqrt(sum(x^2)))
}

# the function from sparcl package
BinarySearch = function(coefs, s) 
{
  if((sum(coefs^2) == 0) | (sum(abs(normalization(coefs))) <= s)) return(0)
  if (s == 1) {return(max(coefs) - 1e-8)}
  lamb1 = 0
  lamb2 = max(abs(coefs)) - 1e-5
  iter = 0
  # while ((iter <= 30) & ((lamb2 - lamb1) > 1e-5)) {
  while ((iter <= 15) & ((lamb2 - lamb1) > 1e-4)) {
    iter = iter + 1
    w_tmp = soft_threshold(coefs, (lamb1 + lamb2) / 2)
    w = normalization(w_tmp)
    if (sum(abs(w)) < s) {
      lamb2 = (lamb1 + lamb2) / 2
    } else {
      lamb1 = (lamb1 + lamb2) / 2
    }
  }
  return((lamb1 + lamb2) / 2)
}

# exact search
ExactSearch = function(coefs, s)
{
  if ((sum(coefs^2) == 0) | (sum(abs(normalization(coefs))) <= s)) return(0)
  if (s == 1) {return(max(coefs) - 1e-8)}
  p = length(coefs)
  ind_seq = (1:p)[1:p != s^2]
  sorted_coefs = sort(coefs, decreasing = TRUE)
  lambda_vec = rep(NA, p)
  for (i in ind_seq) {
    sub_coefs = sorted_coefs[1:i]
    coefs_sum = sum(sub_coefs)
    squared_sum = sum(sub_coefs^2)
    ft = coefs_sum / i
    st_tmp = coefs_sum^2 / i - (coefs_sum^2 - squared_sum * s^2) / (i - s^2)
    if (st_tmp < 0) {warning("The value in square root is negative"); next}
    st = (1 / sqrt(i)) * sqrt(st_tmp)
    # lambda1 = ft + st; lambda2 = ft - st
    # if (((lambda1 + 1e-8) > sorted_coefs[i]) | (lambda1 < (c(sorted_coefs, 0)[i + 1] + 1e-8))) {
    #     lambda1 = NA
    # }
    # if (((lambda2 + 1e-8) > sorted_coefs[i]) | (lambda2 < (c(sorted_coefs, 0)[i + 1] + 1e-8))) {
    #     lambda2 = NA
    # }
    # lambda_vec = c(lambda_vec, c(lambda1, lambda2))
    lambda = ft - st
    if (((lambda + 1e-8) > sorted_coefs[i]) | (lambda < (c(sorted_coefs, 0)[i + 1] + 1e-8))) {
        lambda = NA
    }
    lambda_vec[i] = lambda
  }
  return(min(lambda_vec, na.rm = TRUE))
}




# old version2
# ExactSearch = function(coefs, s)
# {
#   if((sum(coefs^2) == 0) | (sum(abs(normalization(coefs))) <= s)) return(0)
#   p = length(coefs)
#   ind_seq = (1:p)[1:p != s^2]
#   sorted_coefs = sort(coefs, decreasing = TRUE)
#   lambda_vec = rep(NA, p)
#   for (i in ind_seq) {
#     sub_coefs = sorted_coefs[1:i]
#     coefs_sum = sum(sub_coefs)
#     squared_sum = sum(sub_coefs^2)
#     ft = coefs_sum / i
#     st_tmp = coefs_sum^2 / i - (coefs_sum^2 - squared_sum * s^2) / (i - s^2)
#     if (st_tmp < 0) {warning("The value in square root is negative"); next}
#     st = (1 / sqrt(i)) * sqrt(st_tmp)
#     lambda = min(ft + st, ft - st)
#     if ((sorted_coefs[i] > lambda) & (lambda > c(sorted_coefs, 0)[i + 1])) {
#         lambda_vec[i] = lambda
#     }
#   }
#   return(min(lambda_vec, na.rm = TRUE))
# }


# old version
# ExactSearch = function(coefs, s)
# {
#   if((sum(coefs^2) == 0) | (sum(abs(normalization(coefs))) <= s)) return(0)
#   p = length(coefs)
#   ind_seq = (1:p)[1:p != s^2]
#   sorted_coefs = sort(coefs, decreasing = TRUE)
#   lambda_vec = rep(NA, p)
#   for (i in ind_seq) {
#     sub_coefs = sorted_coefs[1:i]
#     coefs_sum = sum(sub_coefs)
#     squared_sum = sum(sub_coefs^2)
#     ft = coefs_sum / i
#     st_tmp = coefs_sum^2 / i - (coefs_sum^2 - squared_sum * s^2) / (i - s^2)
#     if (st_tmp < 0) {warning("st less than 0"); next}
#     st = (1 / sqrt(i)) * sqrt(st_tmp)
#     lambda_vec[i] = min(ft + st, ft - st)
#   }
#   return(min(lambda_vec, na.rm = TRUE))
# }


