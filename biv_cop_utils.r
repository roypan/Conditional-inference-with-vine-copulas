clip_to_unit <- function(value) {
  # clip value to (0, 1).
  pmax(pmin(value, 1 - 1e-7), 1e-7)
}

# Assume starting bivariate copula is permutation summetric.
# u-reflected corresponds to 90 degrees counterclockwise rotation (families 21-30) in VineCopula package
# v-reflected corresponds to 270 degrees counterclockwise rotation (families 31-40) in VineCopula package

# page 272 of Joe (2014): with additional conditional C_{1|2}

# original copula C(u,v)

# survival or doubly-reflected
# Chat_{1|2}(u|v) = 1 - C_{1|2}(1-u|1-v)
# Chat_{2|1}(u|v) = 1 - C_{2|1}(1-v|1-u)

# u-reflected or 1-reflected
# Cacute_{1|2}(u|v) = 1 - C_{1|2}(1-u|v)
# Cacute_{2|1}(v|u) = C_{2|1}(v|1-u)

# v-reflected or 2-reflected
# Cgrave_{1|2}(u|v) = C_{1|2}(u|1-v)
# Cgrave_{2|1}(v|u) = 1 - C_{2|1}(1-v|u)

# VineCopula documentation
# BiCopHfunc1(u1, u2, obj)  # conditional on first argument
# h_1(u_2|u_1,theta) :=  P(U_2 <= u_2 | U_1 = u_1)
# : = partial C(u_1,u_2) / partial u_1, 
# 
# BiCopHfunc2(u1, u2, obj) # conditional on second argument
# h_2(u_1|u_2,theta) := P(U_1 <= u_1 | U_2 = u_2) 
# := partial C(u_1,u_2) / partial u_2, 

# inverse of conditional cdfs 
# BiCopHinv1(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE)
# BiCopHinv2(u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE)
# or
# BiCopHinv1(u1, p , family, par, par2 = 0, obj = NULL, check.pars = TRUE)
# BiCopHinv2(p , u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE)


###########
##  Gaussian-copula  ##
###########

## Gaussian-copula

pbvncop <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=1, par=cpar)
}

logdbvncop <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=1, par=cpar))
}

pcondbvncop <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=1, par=cpar)
}

qcondbvncop <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=1, par=cpar)
}

###########
##  T-copula  ##
###########

## T-copula

pbvtcop <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=2, par=cpar[1], par2=cpar[2])
}

logdbvtcop <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=2, par=cpar[1], par2=cpar[2]))
}

pcondbvtcop <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=2, par=cpar[1], par2=cpar[2])
}

qcondbvtcop <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=2, par=cpar[1], par2=cpar[2])
}


###########
##  MTCJ  ##
###########

## MTCJ

pmtcj <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=3, par=cpar)
}

logdmtcj <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=3, par=cpar))
}

pcond12mtcj <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=3, par=cpar)
}

pcond21mtcj <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=3, par=cpar)
}

qcond12mtcj <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=3, par=cpar)
}

qcond21mtcj <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=3, par=cpar)
}

## MTCJr

pmtcjr <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=13, par=cpar)
}

logdmtcjr <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=13, par=cpar))
}

pcond12mtcjr <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=13, par=cpar)
}

pcond21mtcjr <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=13, par=cpar)
}

qcond12mtcjr <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=13, par=cpar)
}

qcond21mtcjr <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=13, par=cpar)
}

## MTCJu

pmtcju <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=23, par=cpar)
}

logdmtcju <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=23, par=cpar))
}

pcond12mtcju <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=23, par=cpar)
}

pcond21mtcju <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=23, par=cpar)
}

qcond12mtcju <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=23, par=cpar)
}

qcond21mtcju <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=23, par=cpar)
}


## MTCJv

pmtcjv <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=33, par=cpar)
}

logdmtcjv <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=33, par=cpar))
}

pcond12mtcjv <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=33, par=cpar)
}

pcond21mtcjv <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=33, par=cpar)
}

qcond12mtcjv <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=33, par=cpar)
}

qcond21mtcjv <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=33, par=cpar)
}

###########
##  Gumbel  ##
###########

## Gumbel

pgum <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=4, par=cpar)
}

logdgum <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=4, par=cpar))
}

pcond12gum <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=4, par=cpar)
}

pcond21gum <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=4, par=cpar)
}

qcond12gum <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=4, par=cpar)
}

qcond21gum <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=4, par=cpar)
}

## Gumbelr

pgumr <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=14, par=cpar)
}

logdgumr <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=14, par=cpar))
}

pcond12gumr <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=14, par=cpar)
}

pcond21gumr <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=14, par=cpar)
}

qcond12gumr <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=14, par=cpar)
}

qcond21gumr <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=14, par=cpar)
}

## Gumbelu

pgumu <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=24, par=cpar)
}

logdgumu <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=24, par=cpar))
}

pcond12gumu <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=24, par=cpar)
}

pcond21gumu <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=24, par=cpar)
}

qcond12gumu <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=24, par=cpar)
}

qcond21gumu <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=24, par=cpar)
}

## Gumbelv

pgumv <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=34, par=cpar)
}

logdgumv <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=34, par=cpar))
}

pcond12gumv <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=34, par=cpar)
}

pcond21gumv <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=34, par=cpar)
}

qcond12gumv <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=34, par=cpar)
}

qcond21gumv <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=34, par=cpar)
}

###########
##  Frank  ##
###########

## Frank

pfrk <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=5, par=cpar)
}

logdfrk <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=5, par=cpar))
}

pcondfrk <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=5, par=cpar)
}

qcondfrk <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=5, par=cpar)
}

###########
##  Joe  ##
###########

## Joe

pjoe <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=6, par=cpar)
}

logdjoe <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=6, par=cpar))
}

pcond12joe <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=6, par=cpar)
}

pcond21joe <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=6, par=cpar)
}

qcond12joe <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=6, par=cpar)
}

qcond21joe <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=6, par=cpar)
}

## Joer

pjoer <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=16, par=cpar)
}

logdjoer <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=16, par=cpar))
}

pcond12joer <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=16, par=cpar)
}

pcond21joer <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=16, par=cpar)
}

qcond12joer <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=16, par=cpar)
}

qcond21joer <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=16, par=cpar)
}

## Joeu

pjoeu <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=26, par=cpar)
}

logdjoeu <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=26, par=cpar))
}

pcond12joeu <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=26, par=cpar)
}

pcond21joeu <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=26, par=cpar)
}

qcond12joeu <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=26, par=cpar)
}

qcond21joeu <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=26, par=cpar)
}

## Joev

pjoev <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=36, par=cpar)
}

logdjoev <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=36, par=cpar))
}

pcond12joev <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=36, par=cpar)
}

pcond21joev <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=36, par=cpar)
}

qcond12joev <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=36, par=cpar)
}

qcond21joev <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=36, par=cpar)
}


###########
##  BB1  ##
###########

## BB1

pbb1 <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=7, par=cpar[1], par2=cpar[2])
}

logdbb1 <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=7, par=cpar[1], par2=cpar[2]))
}

pcond12bb1 <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=7, par=cpar[1], par2=cpar[2])
}

pcond21bb1 <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=7, par=cpar[1], par2=cpar[2])
}

qcond12bb1 <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=7, par=cpar[1], par2=cpar[2])
}

qcond21bb1 <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=7, par=cpar[1], par2=cpar[2])
}

## BB1r

pbb1r <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=17, par=cpar[1], par2=cpar[2])
}

logdbb1r <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=17, par=cpar[1], par2=cpar[2]))
}

pcond12bb1r <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=17, par=cpar[1], par2=cpar[2])
}

pcond21bb1r <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=17, par=cpar[1], par2=cpar[2])
}

qcond12bb1r <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=17, par=cpar[1], par2=cpar[2])
}

qcond21bb1r <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=17, par=cpar[1], par2=cpar[2])
}

## BB1u

pbb1u <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=27, par=cpar[1], par2=cpar[2])
}

logdbb1u <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=27, par=cpar[1], par2=cpar[2]))
}

pcond12bb1u <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=27, par=cpar[1], par2=cpar[2])
}

pcond21bb1u <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=27, par=cpar[1], par2=cpar[2])
}

qcond12bb1u <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=27, par=cpar[1], par2=cpar[2])
}

qcond21bb1u <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=27, par=cpar[1], par2=cpar[2])
}

## BB1v

pbb1v <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=37, par=cpar[1], par2=cpar[2])
}

logdbb1v <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=37, par=cpar[1], par2=cpar[2]))
}

pcond12bb1v <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=37, par=cpar[1], par2=cpar[2])
}

pcond21bb1v <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=37, par=cpar[1], par2=cpar[2])
}

qcond12bb1v <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=37, par=cpar[1], par2=cpar[2])
}

qcond21bb1v <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=37, par=cpar[1], par2=cpar[2])
}

###########
##  BB8  ##
###########

## BB8

pbb8 <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=10, par=cpar[1], par2=cpar[2])
}

logdbb8 <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=10, par=cpar[1], par2=cpar[2]))
}

pcond12bb8 <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=10, par=cpar[1], par2=cpar[2])
}

pcond21bb8 <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=10, par=cpar[1], par2=cpar[2])
}

qcond12bb8 <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=10, par=cpar[1], par2=cpar[2])
}

qcond21bb8 <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=10, par=cpar[1], par2=cpar[2])
}

## BB8r

pbb8r <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=20, par=cpar[1], par2=cpar[2])
}

logdbb8r <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=20, par=cpar[1], par2=cpar[2]))
}

pcond12bb8r <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=20, par=cpar[1], par2=cpar[2])
}

pcond21bb8r <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=20, par=cpar[1], par2=cpar[2])
}

qcond12bb8r <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=20, par=cpar[1], par2=cpar[2])
}

qcond21bb8r <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=20, par=cpar[1], par2=cpar[2])
}

## BB8u

pbb8u <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=30, par=cpar[1], par2=cpar[2])
}

logdbb8u <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=30, par=cpar[1], par2=cpar[2]))
}

pcond12bb8u <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=30, par=cpar[1], par2=cpar[2])
}

pcond21bb8u <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=30, par=cpar[1], par2=cpar[2])
}

qcond12bb8u <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=30, par=cpar[1], par2=cpar[2])
}

qcond21bb8u <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=30, par=cpar[1], par2=cpar[2])
}

## BB8v

pbb8v <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=40, par=cpar[1], par2=cpar[2])
}

logdbb8v <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=40, par=cpar[1], par2=cpar[2]))
}

pcond12bb8v <- function(u, v, cpar) {
  VineCopula::BiCopHfunc2(u, v, family=40, par=cpar[1], par2=cpar[2])
}

pcond21bb8v <- function(v, u, cpar) {
  VineCopula::BiCopHfunc1(u, v, family=40, par=cpar[1], par2=cpar[2])
}

qcond12bb8v <- function(p, v, cpar) {
  VineCopula::BiCopHinv2(p, v, family=40, par=cpar[1], par2=cpar[2])
}

qcond21bb8v <- function(p, u, cpar) {
  VineCopula::BiCopHinv1(u, p, family=40, par=cpar[1], par2=cpar[2])
}