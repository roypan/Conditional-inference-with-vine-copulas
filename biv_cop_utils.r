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

pmtcju <- function(u, v, cpar) v - pmtcj(1-u, v, cpar)

logdmtcju <- function(u, v, cpar) logdmtcj(1-u, v, cpar)

pcondmtcju <- function(v, u, cpar) pcondmtcj(v, 1-u, cpar)

qcondmtcju <- function(p, u, cpar) qcondmtcj(p, 1-u, cpar)

## MTCJv

pmtcjv <- function(u, v, cpar) u - pmtcj(u, 1-v, cpar)

logdmtcjv <- function(u, v, cpar) logdmtcj(u, 1-v, cpar)

pcondmtcjv <- function(v, u, cpar) 1 - pcondmtcj(1-v, u, cpar)

qcondmtcjv <- function(p, u, cpar) 1 - qcondmtcj(1-p, u, cpar)

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

pgumu <- function(u, v, cpar) v - pgum(1-u, v, cpar)

logdgumu <- function(u, v, cpar) logdgum(1-u, v, cpar)

pcondgumu <- function(v, u, cpar) pcondgum(v, 1-u, cpar)

qcondgumu <- function(p, u, cpar) qcondgum(p, 1-u, cpar)

## Gumbelv

pgumv <- function(u, v, cpar) u - pgum(u, 1-v, cpar)

logdgumv <- function(u, v, cpar) logdgum(u, 1-v, cpar)

pcondgumv <- function(v, u, cpar) 1 - pcondgum(1-v, u, cpar)

qcondgumv <- function(p, u, cpar) 1 - qcondgum(1-p, u, cpar)

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

pjoeu <- function(u, v, cpar) v - pjoe(1-u, v, cpar)

logdjoeu <- function(u, v, cpar) logdjoe(1-u, v, cpar)

pcondjoeu <- function(v, u, cpar) pcondjoe(v, 1-u, cpar)

qcondjoeu <- function(p, u, cpar) qcondjoe(p, 1-u, cpar)

## Joev

pjoev <- function(u, v, cpar) u - pjoe(u, 1-v, cpar)

logdjoev <- function(u, v, cpar) logdjoe(u, 1-v, cpar)

pcondjoev <- function(v, u, cpar) 1 - pcondjoe(1-v, u, cpar)

qcondjoev <- function(p, u, cpar) 1 - qcondjoe(1-p, u, cpar)

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

pbb1u <- function(u, v, cpar) v - pbb1(1-u, v, cpar)

logdbb1u <- function(u, v, cpar) logdbb1(1-u, v, cpar)

pcondbb1u <- function(v, u, cpar) pcondbb1(v, 1-u, cpar)

qcondbb1u <- function(p, u, cpar) qcondbb1(p, 1-u, cpar)

## BB1v

pbb1v <- function(u, v, cpar) u - pbb1(u, 1-v, cpar)

logdbb1v <- function(u, v, cpar) logdbb1(u, 1-v, cpar)

pcondbb1v <- function(v, u, cpar) 1 - pcondbb1(1-v, u, cpar)

qcondbb1v <- function(p, u, cpar) 1 - qcondbb1(1-p, u, cpar)

###########
##  BB6  ##
###########

## BB6

pbb6 <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=8, par=cpar[1], par2=cpar[2])
}

logdbb6 <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=8, par=cpar[1], par2=cpar[2]))
}

pcondbb6 <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=8, par=cpar[1], par2=cpar[2])
}

qcondbb6 <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=8, par=cpar[1], par2=cpar[2])
}

## BB6r

pbb6r <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=18, par=cpar[1], par2=cpar[2])
}

logdbb6r <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=18, par=cpar[1], par2=cpar[2]))
}

pcondbb6r <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=18, par=cpar[1], par2=cpar[2])
}

qcondbb6r <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=18, par=cpar[1], par2=cpar[2])
}

## BB6u

pbb6u <- function(u, v, cpar) v - pbb6(1-u, v, cpar)

logdbb6u <- function(u, v, cpar) logdbb6(1-u, v, cpar)

pcondbb6u <- function(v, u, cpar) pcondbb6(v, 1-u, cpar)

qcondbb6u <- function(p, u, cpar) qcondbb6(p, 1-u, cpar)

## BB6v

pbb6v <- function(u, v, cpar) u - pbb6(u, 1-v, cpar)

logdbb6v <- function(u, v, cpar) logdbb6(u, 1-v, cpar)

pcondbb6v <- function(v, u, cpar) 1 - pcondbb6(1-v, u, cpar)

qcondbb6v <- function(p, u, cpar) 1 - qcondbb6(1-p, u, cpar)

###########
##  BB7  ##
###########

## BB7

pbb7 <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=9, par=cpar[1], par2=cpar[2])
}

logdbb7 <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=9, par=cpar[1], par2=cpar[2]))
}

pcondbb7 <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=9, par=cpar[1], par2=cpar[2])
}

qcondbb7 <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=9, par=cpar[1], par2=cpar[2])
}

## BB7r

pbb7r <- function(u, v, cpar) {
  VineCopula::BiCopCDF(u, v, family=19, par=cpar[1], par2=cpar[2])
}

logdbb7r <- function(u, v, cpar) {
  log(VineCopula::BiCopPDF(u, v, family=19, par=cpar[1], par2=cpar[2]))
}

pcondbb7r <- function(v, u, cpar) {
  VineCopula::BiCopHfunc2(v, u, family=19, par=cpar[1], par2=cpar[2])
}

qcondbb7r <- function(p, u, cpar) {
  VineCopula::BiCopHinv2(p, u, family=19, par=cpar[1], par2=cpar[2])
}

## BB7u

pbb7u <- function(u, v, cpar) v - pbb7(1-u, v, cpar)

logdbb7u <- function(u, v, cpar) logdbb7(1-u, v, cpar)

pcondbb7u <- function(v, u, cpar) pcondbb7(v, 1-u, cpar)

qcondbb7u <- function(p, u, cpar) qcondbb7(p, 1-u, cpar)

## BB7v

pbb7v <- function(u, v, cpar) u - pbb7(u, 1-v, cpar)

logdbb7v <- function(u, v, cpar) logdbb7(u, 1-v, cpar)

pcondbb7v <- function(v, u, cpar) 1 - pcondbb7(1-v, u, cpar)

qcondbb7v <- function(p, u, cpar) 1 - qcondbb7(1-p, u, cpar)

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

pbb8u <- function(u, v, cpar) v - pbb8(1-u, v, cpar)

logdbb8u <- function(u, v, cpar) logdbb8(1-u, v, cpar)

pcondbb8u <- function(v, u, cpar) pcondbb8(v, 1-u, cpar)

qcondbb8u <- function(p, u, cpar) qcondbb8(p, 1-u, cpar)

## BB8v

pbb8v <- function(u, v, cpar) u - pbb8(u, 1-v, cpar)

logdbb8v <- function(u, v, cpar) logdbb8(u, 1-v, cpar)

pcondbb8v <- function(v, u, cpar) 1 - pcondbb8(1-v, u, cpar)

qcondbb8v <- function(p, u, cpar) 1 - qcondbb8(1-p, u, cpar)