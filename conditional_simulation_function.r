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
# Arguments:
# p: vector of length d, e.g. runif(d)
# nsim: sample size for simulation
# A: d*d vine array, or ntrunc*d vine array as only ntrunc rows are used. The variable in the first column is set to the fixed value.
# ntrunc: truncation level between 1 and d-1
# qcond: function for inverse conditional cdf C_{U|V}^{-1}(u|v)
# pcond: function for conditional cdf C_{U|V}(u|v)
# qcondmat: matrix of names of conditional quantile functions for trees 1,...,ntrunc
# pcondmat: matrix of names of conditional cdfs for trees 1,...,ntrunc
# parmat: d*d matrix: for rvinesim1, where all bivariate copula families have 1 parameter, parameter in parmat[ell,j] for ell<j is the parameter of the copula associated with A[ell,j]
# parvec: vector with the union of the parameters associated with the copulas in A[ell,j], j=ell+1,...,d. ell=1,...,ntrunc
# np: d*d matrix of the dimension of the vector for the copulas in A[ell,j], j=ell+1,...,d. ell=1,...,ntrunc; the function will determine parvec[ip1:ip2] for the copula associated with A[ell,j]
# varname: variable name, optional
# extq: the quantile value that the first variable is fixed to
# iprint: print flag for intermediate results
#
# Output: the simulated u-scores with the A[1,1]-th column fixed to extq

rvinesimvec = function(nsim, A, ntrunc, parvec, np, qcondmat, pcondmat, extq, varname = numeric(0), iprint = F) {
  d = ncol(A)
  diagA = diag(A)
  dict = data.frame(Col1=c(0, diagA), Col2=0:d)
  
  A = matrix(dict$Col2[match(A, dict$Col1)], nrow = d, byrow = FALSE)
  
  ii = 0
  ip1 = matrix(0, d, d)
  ip2 = matrix(0, d, d)
  for (ell in 1:ntrunc) {
    for (j in (ell + 1):d) {
      ip1[ell, j] = ii + 1
      ip2[ell, j] = ii + np[ell, j]
      ii = ii + np[ell, j]
    }
  }
  if (iprint) {
    print(ip1)
    print(ip2)
  }
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
  qcond = match.fun(qcondmat[1, 2])
  u[, 2] = qcond(p[, 2], p[, 1], parvec[ip1[1, 2]:ip2[1, 2]])
  qq[, 1, 2] = u[, 2]
  if (icomp[1, 2] == 1) {
    pcond = match.fun(pcondmat[1, 2])
    v[, 1, 2] = pcond(u[, 1], u[, 2], parvec[ip1[1, 2]:ip2[1, 2]])
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
        qcond = match.fun(qcondmat[ell, j])
        qq[, ell, j] = qcond(qq[, ell + 1, j], s, parvec[ip1[ell, j]:ip2[ell, j]])
      }
    }
    qcond = match.fun(qcondmat[1, j])
    qq[, 1, j] = qcond(qq[, 2, j], u[, A[1, j]], parvec[ip1[1, j]:ip2[1, j]])
    u[, j] = qq[, 1, j]
    pcond = match.fun(pcondmat[1, j])
    v[, 1, j] = pcond(u[, A[1, j]], u[, j], parvec[ip1[1, j]:ip2[1, j]])
    if (tt > 1) {
      for (ell in 2:tt) {
        if (A[ell, j] == M[ell, j]) {
          s = qq[, ell, A[ell, j]]
        }
        else {
          s = v[, ell - 1, M[ell, j]]
        }
        if (icomp[ell, j] == 1) {
          pcond = match.fun(pcondmat[ell, j])
          v[, ell, j] = pcond(s, qq[, ell, j], parvec[ip1[ell, j]:ip2[ell, j]])
        }
      }
    }
  }
  
  if (length(varname) != 0) colnames(u) = varname[diagA]
  return(u[, order(dict$Col1[2:(d+1)])])
}
