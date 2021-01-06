library(CopulaModel) # this package in not avaible on CRAN. it can be installed following instructions on https://copula.stat.ubc.ca/
library(VineCopula)

# this function has the same arguments as the rvinesimvec2() function in the CopulaModel R package. please refer to that function for helps and examples
rvinesimvec2_fix1 = function (nsim, A, ntrunc, parvec, np, qcondmat, pcondmat, extq, varname = numeric(0), iprint = F) {
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
	# Note: only the following one line of code is added to the rvinesimvec2 function
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
	return(u)
}
