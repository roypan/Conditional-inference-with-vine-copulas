clip_to_unit <- function(value) {
  # clip value to (0, 1).
  pmax(pmin(value, 1 - 1e-7), 1e-7)
}

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

pcondmtcj <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=3, par=cpar)
}

qcondmtcj <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=3, par=cpar)
}

## MTCJr

pmtcjr <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=13, par=cpar)
}

logdmtcjr <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=13, par=cpar))
}

pcondmtcjr <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=13, par=cpar)
}

qcondmtcjr <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=13, par=cpar)
}

## MTCJu

pmtcju <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=33, par=cpar)
}

logdmtcju <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=33, par=cpar))
}

pcondmtcju <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=33, par=cpar)
}

qcondmtcju <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=33, par=cpar)
}

## MTCJv

pmtcjv <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=23, par=cpar)
}

logdmtcjv <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=23, par=cpar))
}

pcondmtcjv <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=23, par=cpar)
}

qcondmtcjv <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=23, par=cpar)
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

pcondgum <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=4, par=cpar)
}

qcondgum <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=4, par=cpar)
}

## Gumbelr

pgumr <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=14, par=cpar)
}

logdgumr <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=14, par=cpar))
}

pcondgumr <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=14, par=cpar)
}

qcondgumr <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=14, par=cpar)
}

## Gumbelu

pgumu <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=34, par=cpar)
}

logdgumu <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=34, par=cpar))
}

pcondgumu <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=34, par=cpar)
}

qcondgumu <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=34, par=cpar)
}

## Gumbelv

pgumv <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=24, par=cpar)
}

logdgumv <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=24, par=cpar))
}

pcondgumv <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=24, par=cpar)
}

qcondgumv <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=24, par=cpar)
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

pcondjoe <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=6, par=cpar)
}

qcondjoe <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=6, par=cpar)
}

## Joer

pjoer <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=16, par=cpar)
}

logdjoer <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=16, par=cpar))
}

pcondjoer <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=16, par=cpar)
}

qcondjoer <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=16, par=cpar)
}

## Joeu

pjoeu <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=36, par=cpar)
}

logdjoeu <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=36, par=cpar))
}

pcondjoeu <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=36, par=cpar)
}

qcondjoeu <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=36, par=cpar)
}

## Joev

pjoev <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=26, par=cpar)
}

logdjoev <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=26, par=cpar))
}

pcondjoev <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=26, par=cpar)
}

qcondjoev <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=26, par=cpar)
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

pcondbb1 <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=7, par=cpar[1], par2=cpar[2])
}

qcondbb1 <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=7, par=cpar[1], par2=cpar[2])
}

## BB1r

pbb1r <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=17, par=cpar[1], par2=cpar[2])
}

logdbb1r <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=17, par=cpar[1], par2=cpar[2]))
}

pcondbb1r <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=17, par=cpar[1], par2=cpar[2])
}

qcondbb1r <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=17, par=cpar[1], par2=cpar[2])
}

## BB1u

pbb1u <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=37, par=cpar[1], par2=cpar[2])
}

logdbb1u <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=37, par=cpar[1], par2=cpar[2]))
}

pcondbb1u <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=37, par=cpar[1], par2=cpar[2])
}

qcondbb1u <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=37, par=cpar[1], par2=cpar[2])
}

## BB1v

pbb1v <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=27, par=cpar[1], par2=cpar[2])
}

logdbb1v <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=27, par=cpar[1], par2=cpar[2]))
}

pcondbb1v <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=27, par=cpar[1], par2=cpar[2])
}

qcondbb1v <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=27, par=cpar[1], par2=cpar[2])
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

pcondbb8 <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=10, par=cpar[1], par2=cpar[2])
}

qcondbb8 <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=10, par=cpar[1], par2=cpar[2])
}

## BB8r

pbb8r <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=20, par=cpar[1], par2=cpar[2])
}

logdbb8r <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=20, par=cpar[1], par2=cpar[2]))
}

pcondbb8r <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=20, par=cpar[1], par2=cpar[2])
}

qcondbb8r <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=20, par=cpar[1], par2=cpar[2])
}

## BB8u

pbb8u <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=40, par=cpar[1], par2=cpar[2])
}

logdbb8u <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=40, par=cpar[1], par2=cpar[2]))
}

pcondbb8u <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=40, par=cpar[1], par2=cpar[2])
}

qcondbb8u <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=40, par=cpar[1], par2=cpar[2])
}

## BB8v

pbb8v <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=30, par=cpar[1], par2=cpar[2])
}

logdbb8v <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=30, par=cpar[1], par2=cpar[2]))
}

pcondbb8v <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=30, par=cpar[1], par2=cpar[2])
}

qcondbb8v <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=30, par=cpar[1], par2=cpar[2])
}