library(VineCopula)
source('biv_cop_utils.r')
source('conditional_simulation_function.r')

set.seed(10)

d = 4 # dimension
param1 = matrix(c(0, .8,.6,.7, 0, 0, 1.5,1.5, 0, 0, 0, .3, 0, 0, 0, 0) , d, d, byrow = TRUE) # parameter 1 array
#     [,1] [,2] [,3] [,4]
#[1,]    0  0.8  0.6  0.7
#[2,]    0  0.0  1.5  1.5
#[3,]    0  0.0  0.0  0.3
#[4,]    0  0.0  0.0  0.0

param2 = matrix(0, d, d) # parameter 2 array

# variable 1 is set to a fixed value
D = matrix(c(1,1,2,3,0,2,1,2,0,0,3,1,0,0,0,4), d, d, byrow = TRUE) # vine array
#     [,1] [,2] [,3] [,4]
#[1,]    1    1    2    3
#[2,]    0    2    1    2
#[3,]    0    0    3    1
#[4,]    0    0    0    4

fam = matrix(c(0, 1, 1, 1, 0, 0, 5, 5, 0, 0, 0, 5, 0, 0, 0, 0), d, d, byrow = TRUE) # family array
#     [,1] [,2] [,3] [,4]
#[1,]    0    1    1    1
#[2,]    0    0    5    5
#[3,]    0    0    0    5
#[4,]    0    0    0    0


rvinesimvec(5, D, ntrunc = 3, fam, param1, param2, extq = .95, iprint = FALSE)
#    [,1]      [,2]      [,3]      [,4]
#[1,] 0.95 0.8060706 0.8696160 0.7819850
#[2,] 0.95 0.8305604 0.8432464 0.3824218
#[3,] 0.95 0.8295451 0.4648344 0.4259508
#[4,] 0.95 0.9322303 0.8989036 0.8293215
#[5,] 0.95 0.8867755 0.7618500 0.9198797



set.seed(20)
# variable 2 is set to a fixed value; the vine structure is the same as the first case
D = matrix(c(2,2,2,3,0,1,1,2,0,0,3,1,0,0,0,4), d, d, byrow = TRUE) # vine array
#     [,1] [,2] [,3] [,4]
#[1,]    2    2    2    3
#[2,]    0    1    1    2
#[3,]    0    0    3    1
#[4,]    0    0    0    4
param1 = matrix(c(0, .8,.6,.7, 0, 0, 1.5,1.5, 0, 0, 0, .3, 0, 0, 0, 0) , d, d, byrow = TRUE) # parameter 1 array

rvinesimvec(5, D, ntrunc = 3, fam, param1, param2, extq = .95, iprint = FALSE)
#          [,1] [,2]      [,3]       [,4]
#[1,] 0.9946533 0.95 0.9598202 0.91564331
#[2,] 0.6971869 0.95 0.8922886 0.76801945
#[3,] 0.6677627 0.95 0.0721232 0.05152934
#[4,] 0.8526675 0.95 0.9143433 0.77864483
#[5,] 0.8679792 0.95 0.5914071 0.85432253



set.seed(30)
# variable 3 is set to a fixed value; the vine structure is the same as the first case
D = matrix(c(3,3,3,2,0,2,2,3,0,0,4,4,0,0,0,1), d, d, byrow = TRUE) # vine array
#     [,1] [,2] [,3] [,4]
#[1,]    3    3    3    2
#[2,]    0    2    2    3
#[3,]    0    0    4    4
#[4,]    0    0    0    1
param1 = matrix(c(0, .6,.7,.8, 0, 0, 1.5,1.5, 0, 0, 0, .3, 0, 0, 0, 0) , d, d, byrow = TRUE) # parameter 1 array

rvinesimvec(5, D, ntrunc = 3, fam, param1, param2, extq = .95, iprint = FALSE)
#[1,] 0.8660866 0.5594631 0.95 0.4699688
#[2,] 0.9061612 0.9775554 0.95 0.8904169
#[3,] 0.8595563 0.6475531 0.95 0.8530475
#[4,] 0.9861474 0.9927880 0.95 0.9865439
#[5,] 0.5902382 0.5503248 0.95 0.6616354



set.seed(40)
# variable 4 is set to a fixed value; the vine structure is the same as the first case
D = matrix(c(4,4,3,2,0,3,4,3,0,0,2,4,0,0,0,1), d, d, byrow = TRUE) # vine array
#     [,1] [,2] [,3] [,4]
#[1,]    4    4    3    2
#[2,]    0    3    4    3
#[3,]    0    0    2    4
#[4,]    0    0    0    1
param1 = matrix(c(0, .7,.6,.8, 0, 0, 1.5,1.5, 0, 0, 0, .3, 0, 0, 0, 0) , d, d, byrow = TRUE) # parameter 1 array

rvinesimvec(5, D, ntrunc = 3, fam, param1, param2, extq = .95, iprint = FALSE)
#[1,] 0.7357954 0.4537622 0.8603368 0.95
#[2,] 0.7605298 0.8548272 0.7120318 0.95
#[3,] 0.8985691 0.4572681 0.9058219 0.95
#[4,] 0.3688029 0.6139390 0.8216701 0.95
#[5,] 0.9289678 0.9670506 0.6495283 0.95

# Remark: The general algorithm to convert to a vine array with a given variable in the first column will be presented elsewhere