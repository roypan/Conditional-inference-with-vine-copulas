# the function to perform cross prediction for the index-th variable
# udat: a n*d matrix of u-scores, the index-th column can be NA or any arbitrary values
# RVM: an RVM object in the VineCopula package
# index: the column index of the response variable
# q: the quantile to be predicted; 0.5 leads to median
#
# Output: returns the predicted u-cores for the index-th variable

RVineCrossPrediction = function(udat, RVM, index, q) {
  sol = c()
  n_var = ncol(udat)
  n_obs = nrow(udat)
  
  gl_object = gauss.quad.prob(30, dist = 'uniform')
  
  for (i in 1:n_obs) {
    joint_density = function(u) { #must accept vector input
      dat = t(replicate(length(u), udat[i, ]))
      dat[, index] = u
      l = RVineLogLik(dat, RVM, separate = TRUE)
      return(exp(l$loglik))
    }
    
    dist_denum = sum(gl_object$weights * joint_density(gl_object$nodes))
    
    obj = function(x) {
      return(x * sum(gl_object$weights * joint_density(x * gl_object$nodes)) - dist_denum * q) # compare with integrate() # use try function: <<-
    }
    sol = c(sol, uniroot(obj, c(0, 1))$root)
  }
  
  return(sol)
}
