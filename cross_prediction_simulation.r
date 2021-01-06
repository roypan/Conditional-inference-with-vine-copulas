source('cross_prediction.r')

N = 100
n_train = 800
n_test = 200
n_var = 5
alpha = .2


# example of using the cross prediction function
# define 5-dimensional R-vine tree structure matrix
Matrix = c(5, 1, 4, 2, 3,
            0, 4, 1, 3, 2,
            0, 0, 3, 1, 2,
            0, 0, 0, 2, 1,
            0, 0, 0, 0, 1)
Matrix = matrix(Matrix, 5, 5)

# define R-vine pair-copula family matrix
family = c(0, 1, 1, 1, 1,
            0, 0, 1, 1, 1,
            0, 0, 0, 1, 1,
            0, 0, 0, 0, 1,
            0, 0, 0, 0, 0)
family = matrix(family, 5, 5)

# define R-vine pair-copula parameter matrix
par = c(0, .1, .2, .3, .7,
         0, 0, .3, .5, .5,
         0, 0, 0, .4, .6,
         0, 0, 0, 0, .8,
         0, 0, 0, 0, 0)
par = matrix(par, 5, 5)

par2 = matrix(0, 5, 5)

RVM = RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3", "V4", "V5"))

set.seed(10)

for (sim_num in 1:N) {
	print(sim_num)
	simdata_u = RVineSim(1000, RVM)
	simdata_x = qnorm(simdata_u)
	
	simdata_u_train = simdata_u[1:n_train, ]
	simdata_u_test = simdata_u[-(1:n_train), ]
	simdata_x_train = simdata_x[1:n_train, ]
	simdata_x_test = simdata_x[-(1:n_train), ]
	
	# gaussian copula
	RVM = RVineStructureSelect(simdata_u_train, familyset  = 1, progress = TRUE) 
	
	pred_results_gc_median = matrix(nrow = n_test, ncol = n_var)
	pred_results_gc_lb = matrix(nrow = n_test, ncol = n_var)
	pred_results_gc_ub = matrix(nrow = n_test, ncol = n_var)
	
	for (i in 1:n_var) {
		pred_results_gc_median[, i] = qnorm(RVineCrossPrediction(simdata_u_test, RVM, i, .5))
		pred_results_gc_lb[, i] = qnorm(RVineCrossPrediction(simdata_u_test, RVM, i, .1))
		pred_results_gc_ub[, i] = qnorm(RVineCrossPrediction(simdata_u_test, RVM, i, .9))
	}
	
	# vine copula
	RVM = RVineStructureSelect(simdata_u_train, familyset  = 1:10, progress = TRUE) 
	
	pred_results_vc_median = matrix(nrow = n_test, ncol = n_var)
	pred_results_vc_lb = matrix(nrow = n_test, ncol = n_var)
	pred_results_vc_ub = matrix(nrow = n_test, ncol = n_var)
	
	for (i in 1:n_var) {
		pred_results_vc_median[, i] = qnorm(RVineCrossPrediction(simdata_u_test, RVM, i, .5))
		pred_results_vc_lb[, i] = qnorm(RVineCrossPrediction(simdata_u_test, RVM, i, .1))
		pred_results_vc_ub[, i] = qnorm(RVineCrossPrediction(simdata_u_test, RVM, i, .9))
	}
}
