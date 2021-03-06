# function to convert a vine arrat A to a maximum matrix
# This is Algorithm 5 in Joe (2014)
# varray2M : maximum array computed from vine array A, used for log-likelihood and simulation
# A = dxd vine array of R-vine in standard order with 1:d on diagonal
# iprint = print flag for intermediate calculations
# str = string to describe vine for printing if iprint=T
# Output: list with two components  
#   M = dxd array with m_{kj}= max a_{k1},..,a_{kj}
#   icomp = dxd indicator array on whether back step [k,j] is needed
#   icomp[k-1,m_{kj}=1 if  a_{kj}<m_{kj} for k>=2

varray2M = function(A, iprint = F, str = "") {
  d = ncol(A)
  d1 = d - 1
  M = A
  icomp = matrix(0, d, d)
  for (k in 2:d1) {
    for (j in (k + 1):d) M[k, j] = max(M[k - 1, j], A[k, 
      j])
  }
  if (iprint) {
    cat("\n", str, "\n")
    print(A)
    print(M)
  }
  for (k in 2:d1) {
    for (j in (k + 1):d) {
      if (A[k, j] < M[k, j]) 
        icomp[k - 1, M[k, j]] = 1
    }
  }
  if (iprint) 
    print(icomp)
  list(mxarray = M, icomp = icomp)
}


################################################################################################################################################
# simulation for R-vine copulas when the first variable is set to a fixed value
#
# bivariate copula C_{12}(u,v) = pcop(u,v) with conditional distributions
#   C_{2|1}(v|u)=pcond21(v,u), inverse C^{-1}_{2|1}(p|u)= qcond21(p,u)
#   C_{1|2}(v|u)=pcond12(v,u), inverse C^{-1}_{1|2}(p|v)= qcond21(p,v)
#
# Arguments:
# p: vector of length d, e.g. runif(d)
# nsim: sample size for simulation
# A: d*d vine array, or ntrunc*d vine array as only ntrunc rows are used. The variable in the first column is set to the fixed value.
# ntrunc: truncation level between 1 and d-1
# fam: d by d matrix for the copula families coded by the VineCopula package
# param1: d by d matrix for the first parameter of the copula models
# param2: d by d matrix for the second parameter of the copula models
# varname: variable name, optional
# extq: the quantile value that the first variable is fixed to
# iprint: print flag for intermediate results
#
# Output: the simulated u-scores with the A[1,1]-th column fixed to extq

rvinesimvec = function(nsim, A, ntrunc, fam, param1, param2, extq, varname = numeric(0), iprint = F) {
  d = ncol(A)
  diagA = diag(A)
  dict = data.frame(Col1=c(0, diagA), Col2=0:d)
  
  # The next line temporarily permutes variable indices so that diagonal of vine array is 1:d in order to apply Algorithm 17 in Joe (2014)
  A = matrix(dict$Col2[match(A, dict$Col1)], nrow = d, byrow = FALSE)
  
  out = varray2M(A)
  M = out$mxarray
  icomp = out$icomp
  p = matrix(runif(nsim * d), nsim, d)
  # Note: only the following one line of code is modified compared with Algorithm 17 in Joe (2014)
  p[, 1] = rep(extq, nsim)
  qq = array(0, c(nsim, d, d))
  v = array(0, c(nsim, d, d))
  u = matrix(0, nsim, d)
  u[, 1] = p[, 1]
  qq[, 1, 1] = p[, 1]
  qq[, 2, 2] = p[, 2]
  u[, 2] = BiCopHinv1(p[, 1], p[, 2], family = fam[1, 2], par = param1[1, 2], par2 = param2[1, 2])
  qq[, 1, 2] = u[, 2]
  if (icomp[1, 2] == 1) {
    v[, 1, 2] = BiCopHfunc2(u[, 1], u[, 2], family = fam[1, 2], par = param1[1, 2], par2 = param2[1, 2])
  }
  for (j in 3:d) {
    tt = min(ntrunc, j - 1)
    qq[, tt + 1, j] = p[, j]
    if (tt > 1) {
      for (ell in seq(tt, 2)) {
        if (A[ell, j] == M[ell, j]) {
          s = qq[, ell, A[ell, j]]
        }
        else {
          s = v[, ell - 1, M[ell, j]]
        }
        qq[, ell, j] = BiCopHinv1(s, qq[, ell + 1, j], family = fam[ell, j], par = param1[ell, j], par2 = param2[ell, j])
      }
    }
    qq[, 1, j] = BiCopHinv1(u[, A[1, j]], qq[, 2, j], family = fam[1, j], par = param1[1, j], par2 = param2[1, j])
    u[, j] = qq[, 1, j]
    v[, 1, j] = BiCopHfunc2(u[, A[1, j]], u[, j], family = fam[1, j], par = param1[1, j], par2 = param2[1, j])
    if (tt > 1) {
      for (ell in 2:tt) {
        if (A[ell, j] == M[ell, j]) {
          s = qq[, ell, A[ell, j]]
        }
        else {
          s = v[, ell - 1, M[ell, j]]
        }
        if (icomp[ell, j] == 1) {
          v[, ell, j] = BiCopHfunc2(s, qq[, ell, j], family = fam[ell, j], par = param1[ell, j], par2 = param2[ell, j])
        }
      }
    }
  }
  
  if (length(varname) != 0) colnames(u) = varname[diagA]
  return(u[, order(dict$Col1[2:(d+1)])])
}
