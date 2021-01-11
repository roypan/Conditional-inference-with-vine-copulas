source('conditional_simulation.r')

d=4
np=matrix(0,d,d)
np[1,2:d]=1; np[2,3:d]=1; np[3,4:d]=1
qcondnames=c("qcondgum","qcondfrk","qcondfrk")
pcondnames=c("pcondgum","pcondfrk","pcondfrk")
parvec=c(.8,.6,.7, 1.5,1.5, .3)
D = matrix(c(1,1,2,3,0,2,1,2,0,0,3,1,0,0,0,4), d, d, byrow = TRUE)
pcondmat = matrix(c("",rep("pcondbvncop",3),"","",rep("pcondfrk",2),"","","","pcondfrk",rep("",4)), 4, 4, byrow=TRUE)
qcondmat = matrix(c("",rep("qcondbvncop",3),"","",rep("qcondfrk",2),"","","","qcondfrk",rep("",4)), 4, 4, byrow=TRUE)

rvinesimvec2_fix1(100, D, ntrunc = 3, parvec, np, qcondmat, pcondmat, extq = .95, iprint = FALSE)
