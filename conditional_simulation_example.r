library(VineCopula)
source('biv_cop_utils.r')
source('conditional_simulation_function.r')

set.seed(10)

d = 4 # dimension
np = matrix(0,d,d) # number of parameters on each edge
np[1,2:d] = 1; np[2,3:d] = 1; np[3,4:d] = 1 
parvec1 = c(.8,.6,.7, 1.5,1.5, .3) # parameter vector

# variable 1 is set to a fixed value
D = matrix(c(1,1,2,3,0,2,1,2,0,0,3,1,0,0,0,4), d, d, byrow = TRUE) # vine array
#     [,1] [,2] [,3] [,4]
#[1,]    1    1    2    3
#[2,]    0    2    1    2
#[3,]    0    0    3    1
#[4,]    0    0    0    4

pcondmat = matrix(c("",rep("pcondbvncop",3),"","",rep("pcondfrk",2),"","","","pcondfrk",rep("",4)), 4, 4, byrow=TRUE)
#     [,1] [,2]          [,3]          [,4]         
#[1,] ""   "pcondbvncop" "pcondbvncop" "pcondbvncop"
#[2,] ""   ""            "pcondfrk"    "pcondfrk"   
#[3,] ""   ""            ""            "pcondfrk"   
#[4,] ""   ""            ""            ""        
qcondmat = matrix(c("",rep("qcondbvncop",3),"","",rep("qcondfrk",2),"","","","qcondfrk",rep("",4)), 4, 4, byrow=TRUE)
#     [,1] [,2]          [,3]          [,4]         
#[1,] ""   "qcondbvncop" "qcondbvncop" "qcondbvncop"
#[2,] ""   ""            "qcondfrk"    "qcondfrk"   
#[3,] ""   ""            ""            "qcondfrk"   
#[4,] ""   ""            ""            ""           


rvinesimvec(5, D, ntrunc = 3, parvec1, np, qcondmat, pcondmat, extq = .95, iprint = FALSE)
#    [,1]      [,2]      [,3]      [,4]
#[1,] 0.95 0.8060706 0.8696160 0.7819850
#[2,] 0.95 0.8305604 0.8432464 0.3824218
#[3,] 0.95 0.8295451 0.4648344 0.4259508
#[4,] 0.95 0.9322303 0.8989036 0.8293215
#[5,] 0.95 0.8867755 0.7618500 0.9198797


set.seed(10)
D = matrix(c(4,4,1,3,0,1,4,1,0,0,3,4,0,0,0,2), d, d, byrow = TRUE) # the same vine array as above, except indices (1,2,3,4) are permuted to (4,1,3,2)
#     [,1] [,2] [,3] [,4]
#[1,]    4    4    1    3
#[2,]    0    1    4    1
#[3,]    0    0    3    4
#[4,]    0    0    0    2
rvinesimvec(5, D, ntrunc = 3, parvec1, np, qcondmat, pcondmat, extq = .95, iprint = FALSE) # generated data are the same
#          [,1]      [,2]      [,3] [,4]
#[1,] 0.8060706 0.7819850 0.8696160 0.95
#[2,] 0.8305604 0.3824218 0.8432464 0.95
#[3,] 0.8295451 0.4259508 0.4648344 0.95
#[4,] 0.9322303 0.8293215 0.8989036 0.95
#[5,] 0.8867755 0.9198797 0.7618500 0.95



set.seed(20)
# variable 2 is set to a fixed value; the vine structure is the same as the first case
D = matrix(c(2,2,2,3,0,1,1,2,0,0,3,1,0,0,0,4), d, d, byrow = TRUE) # vine array
#     [,1] [,2] [,3] [,4]
#[1,]    2    2    2    3
#[2,]    0    1    1    2
#[3,]    0    0    3    1
#[4,]    0    0    0    4
parvec2 = c(.8,.6,.7, 1.5,1.5, .3) # parameter vector

rvinesimvec(5, D, ntrunc = 3, parvec2, np, qcondmat, pcondmat, extq = .95, iprint = FALSE)
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
parvec3 = c(.6,.7,.8, 1.5,1.5, .3) # parameter vector

rvinesimvec(5, D, ntrunc = 3, parvec3, np, qcondmat, pcondmat, extq = .95, iprint = FALSE)
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
parvec4 = c(.7,.6,.8, 1.5,1.5, .3) # parameter vector

rvinesimvec(5, D, ntrunc = 3, parvec4, np, qcondmat, pcondmat, extq = .95, iprint = FALSE)
#[1,] 0.7357954 0.4537622 0.8603368 0.95
#[2,] 0.7605298 0.8548272 0.7120318 0.95
#[3,] 0.8985691 0.4572681 0.9058219 0.95
#[4,] 0.3688029 0.6139390 0.8216701 0.95
#[5,] 0.9289678 0.9670506 0.6495283 0.95