library(VineCopula)

################################################################################################################################################
# utility functions

# function to convert a vine arrat A to a maximum  =matrix
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


# bivariate copula conditional quantile functions
qcondbb2 = function(p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) {
    th = cpar[1]
    de = cpar[2]
    de1 = 1/de
    th1 = 1/th
    ut = u^(-th) - 1
    x = exp(ut * de) - 1
    den = (1 + ut)^(-th1 - 1)/(x + 1)
    if (is.infinite(x) | den <= 1e-100) {
        if (iprint) 
            cat("\n** infinite x ", p, u, th, de, ut, "\n")
        lx = ut * de
        r = 1
        tha = -th/(1 + th)
        diff = 1
        iter = 0
        while (abs(diff) > eps & iter < mxiter) {
            sm = 1 + r
            smd = de1 * log(sm)
            rhs = (p * sm)^tha
            g = 1 + log(sm)/lx - rhs
            gp = 1/sm/lx - tha * rhs/sm
            iter = iter + 1
            diff = g/gp
            r = r - diff
            while (r <= 0) {
                diff = diff/2
                r = r + diff
            }
            if (iprint) 
                cat(iter, r, diff, "\n")
        }
        v = (1 + de1 * (log(r) + lx))^(-th1)
        if (iprint) 
            cat("** infinite x ", p, u, th, de, ut, v, 
                "\n")
        return(v)
    }
    pden = p * den
    rhs = (1 + log(2 * x + 1))^(-1 - th1)/pden
    y = rhs - 1 + x
    if (y <= 0) 
        y = 0.1
    if (y > 1000) 
        y = 1000
    diff = 1
    iter = 0
    while ((abs(diff/y) > eps) & iter < mxiter) {
        sm = x + y + 1
        smd = de1 * log(sm)
        G21 = ((1 + smd)^(-th1 - 1))/sm
        gpdf = -G21
        gpdf = gpdf/(1 + smd)/sm/de/th
        gpdf = gpdf * (1 + th + th * de * (1 + smd))
        iter = iter + 1
        diff = (G21 - pden)/gpdf
        y = y - diff
        while (y <= 0) {
            diff = diff/2
            y = y + diff
        }
        if (iprint) 
            cat(iter, y, diff, "\n")
    }
    v = (1 + de1 * log(y + 1))^(-th1)
    if (iprint) 
        cat(v, "ended at iter. ", iter, "\n")
    v
}

qcondbb3 = function(p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) {
    th = cpar[1]
    de = cpar[2]
    de1 = 1/de
    th1 = 1/th
    dt = (de^th1)
    ul = -log(u)
    ut = (ul^th)
    x = exp(ut * de) - 1
    if (is.infinite(x) | x > 1e+200) {
        if (iprint) 
            cat("\n** infinite x ", p, u, th, de, ut, "\n")
        lx = ut * de
        llx = log(lx)
        lxt = (lx^th1)
        con = log(p) - lx + (th1 - 1) * llx - lxt/dt
        r = 1
        diff = 1
        iter = 0
        while (abs(diff) > eps & iter < mxiter) {
            sm = r + 1
            sml = log(sm) + lx
            smlt = (sml^th1)
            h = sml + (1 - th1) * log(sml) + smlt/dt + con
            hp = 1/sm + (1 - th1)/sml/sm + th1 * smlt/sml/sm/dt
            iter = iter + 1
            diff = h/hp
            r = r - diff
            while (r <= 0) {
                diff = diff/2
                r = r + diff
            }
            if (iprint) 
                cat(iter, r, diff, "\n")
        }
        v = ((log(r) + lx)^th1)/dt
        v = exp(-v)
        if (iprint) 
            cat("** infinite x ", p, u, th, de, ut, v, 
                "\n")
        return(v)
    }
    if (x <= 1e-10) {
        if (iprint) 
            cat("\n** x near 0 ", p, u, th, de, ut, "\n")
        xx = ut * de
        xxt = (xx^th1)
        con = log(p) - xxt/dt
        r = 1
        diff = 1
        iter = 0
        while (abs(diff) > eps & iter < mxiter) {
            sm = r + 1
            sml = log(sm)
            smt = (sm^th1)
            h = xx * r + (1 - th1) * sml + smt * xxt/dt + con
            hp = xx + (1 - th1)/sm + th1 * smt * xxt/sm/dt
            iter = iter + 1
            diff = h/hp
            r = r - diff
            while (r <= 0) {
                diff = diff/2
                r = r + diff
            }
            if (iprint) 
                cat(iter, r, diff, "\n")
        }
        v = (r * xx/de)^th1
        v = exp(-v)
        if (iprint) 
            cat("** x near 0 ", p, u, th, de, ut, v, "\n")
        return(v)
    }
    lx = ut * de
    llx = log(lx)
    lxt = (lx^th1)
    con = log(p) - lx + (th1 - 1) * llx - lxt/dt
    y = x
    if (y <= 0) 
        y = 0.1
    if (y > 1000) 
        y = 1000
    if (iprint) 
        cat("p, u ", p, u, "\n")
    if (iprint) 
        cat("x con and starting y: ", x, con, y, "\n")
    diff = 1
    iter = 0
    while ((abs(diff/y) > eps) & iter < mxiter) {
        sm = x + y + 1
        sml = log(sm)
        smlt = (sml^th1)
        h = sml + (1 - th1) * log(sml) + smlt/dt + con
        hp = 1/sm + (1 - th1)/sml/sm + th1 * smlt/sml/sm/dt
        iter = iter + 1
        diff = h/hp
        y = y - diff
        while (y <= 0) {
            diff = diff/2
            y = y + diff
        }
        if (iprint) 
            cat(iter, y, diff, "\n")
    }
    v = (de1 * log(y + 1))^th1
    v = exp(-v)
    if (iter >= mxiter) {
        cat("** did not converge **\n")
        cat("p,u,th,de,lasty,v: ", p, u, th, de, y, v, 
            "\n")
    }
    if (iprint) 
        cat("v and #iter :", v, iter, "\n")
    v
}

qcondbb4 = function(p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) {
    th = cpar[1]
    de = cpar[2]
    de1 = 1/de
    th1 = 1/th
    be = 4 * pbb4(0.5, 0.5, cpar) - 1
    ut = u^(-th) - 1
    x = ut^(-de)
    con = (de1 + 1) * log(x) - (th + 1) * log(u) - log(p)
    iter = 0
    diff = 1
    y = 0.5 * x
    if (de > 1.7 | u > 0.9 | p > 0.9) {
        diff = 1
        v = u
        v1 = 0
        v2 = 1
        vold = v
        if (iprint) 
            cat("\n (p,u)=", p, u, "\n")
        while (diff > eps) {
            ccdf = pcondbb4(v, u, cpar)
            di = ccdf - p
            if (di < 0) {
                v1 = v
            }
            else {
                v2 = v
            }
            diff = v2 - v1
            vold = v
            v = (v1 + v2)/2
            if (iprint) 
                cat(vold, ccdf, di, "\n")
        }
        return(v)
    }
    if (th < 0.2) {
        v = qcondgal(p, u, de)
        vt = v^(-th) - 1
        y = vt^(-de)
    }
    if (be >= 0.8) 
        y = x
    while (iter < mxiter & max(abs(diff)) > eps) {
        xy = x + y
        vt = y^(-de1)
        tem = xy^(-de1)
        sm = ut + vt + 1 - tem
        xtem = ut/x - tem/xy
        ytem = vt/y - tem/xy
        h = -(th1 + 1) * log(sm) + log(xtem) + con
        hp = (th1 + 1) * ytem/sm/de + (1 + de1) * tem/xtem/xy/xy
        diff = h/hp
        y = y - diff
        if (iprint) 
            cat(iter, diff, y, "\n")
        while (min(y) <= 0 | max(abs(diff)) > 5) {
            diff = diff/2
            y = y + diff
        }
        iter = iter + 1
    }
    if (iprint & iter >= mxiter) 
        cat("***did not converge\n")
    (1 + y^(-de1))^(-th1)
}

qcondbb5 = function(p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) {
    th = cpar[1]
    de = cpar[2]
    de1 = 1/de
    th1 = 1/th
    x = -log(u)
    xt = x^th
    xdt = x^(-de * th)
    con = x - log(p)
    iter = 0
    diff = 1
    y = 0.5 * x
    be = 4 * pbb5(0.5, 0.5, cpar) - 1
    if (be >= 0.8) 
        y = x
    while (iter < mxiter & max(abs(diff)) > eps) {
        yt = y^th
        ydt = y^(-de * th)
        xydt = xdt + ydt
        xyp = xydt^(-1/de)
        w = xt + yt - xyp
        wth = w^(1/th)
        zx = (xt - xyp/xydt * xdt)/x
        zy = (yt - xyp/xydt * ydt)/y
        h = -wth + (th1 - 1) * log(w) + log(zx) + con
        hp = (-wth + (1 - th)) * zy/w - th * (1 + de) * (xyp/xydt/xydt) * 
            (xdt/x) * (ydt/y)/zx
        diff = h/hp
        y = y - diff
        if (iprint) 
            cat(iter, diff, y, "\n")
        while (min(y) <= 0 | max(abs(diff)) > 5) {
            diff = diff/2
            y = y + diff
        }
        iter = iter + 1
    }
    if (iprint & iter >= mxiter) 
        cat("***did not converge\n")
    exp(-y)
}

qcondbb6 = function(p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) {
    th = cpar[1]
    de = cpar[2]
    de1 = 1/de
    th1 = 1/th
    be = 4 * pbb6(0.5, 0.5, cpar) - 1
    ubar = 1 - u
    zu = 1 - ubar^th
    x = -log(zu)
    xd = x^de
    con = (de - 1) * log(x) - (th1 - 1) * log(1 - zu) + x - log(p)
    iter = 0
    diff = 1
    y = 0.5 * x
    if (be >= 0.8) 
        y = x
    while (iter < mxiter & max(abs(diff)) > eps) {
        yd = y^de
        sm = xd + yd
        tem = sm^(1/de)
        w = exp(-tem)
        h = (th1 - 1) * log(1 - w) - tem + (de1 - 1) * log(sm) + 
            con
        hp = (th1 - 1) * w * tem * yd/(1 - w)/sm/y - tem * yd/sm/y + 
            (1 - de) * yd/y/sm
        diff = h/hp
        y = y - diff
        if (iprint) 
            cat(iter, diff, y, "\n")
        if (any(is.nan(y))) {
            if (iprint) 
                cat(p, u, cpar, x, "\n")
            return(u)
        }
        while (min(y) <= 0 | max(abs(diff)) > 5) {
            diff = diff/2
            y = y + diff
        }
        iter = iter + 1
    }
    if (iprint & iter >= mxiter) 
        cat("***did not converge\n")
    1 - (1 - exp(-y))^th1
}

qcondbb7 = function(p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) {
    th = cpar[1]
    de = cpar[2]
    de1 = 1/de
    th1 = 1/th
    ut = 1 - (1 - u)^th
    x = ut^(-de) - 1
    if (x <= 0) {
        v = 0.999999
        if (u > v) 
            v = u
        if (iprint) 
            cat("\n**** x below 1 ", p, u, x, "\n")
        return(v)
    }
    den = (1 - ut)^(th1 - 1) * ut/(x + 1)
    pden = p * den
    rhs = pden * (1 - (2 * x + 1)^(-de1))^(1 - th1)
    y = rhs^(-de/(de + 1)) - 1 - x
    if (y <= 0) 
        y = 0.1
    if (th > 3 & (u > 0.8 | p > 0.9)) {
        diff = 1
        v = u
        v1 = 0
        v2 = 1
        vold = v
        if (iprint) 
            cat("\n (p,u)=", p, u, "\n")
        while (diff > eps) {
            ccdf = pcondbb7(v, u, cpar)
            di = ccdf - p
            if (di < 0) {
                v1 = v
            }
            else {
                v2 = v
            }
            diff = v2 - v1
            vold = v
            v = (v1 + v2)/2
            if (iprint) 
                cat(vold, ccdf, di, "\n")
        }
        return(v)
    }
    if (x < 1e-05) {
        epsx = de * (1 - ut)
        tem = p * (1 - (1 + de1) * epsx)
        tem = tem^(-th/(th - 1)) - 1
        epsy = tem * epsx
        if (epsy > 1e-05) 
            epsy = 1e-05
        y = epsy
        if (iprint) 
            cat("\nsmall x case ", x, y, "\n")
    }
    if (th < 1.01) {
        thr = -th/(1 + th)
        tem = (p^thr) - 1
        tem = tem * (u^(-th)) + 1
        y = (tem - 1)^de
    }
    else if (de < 0.1) {
        v = 1 - qcondjoe(p, u, de)
        y = v^th
        y = (1 - y)^(-de) - 1
    }
    diff = 1
    iter = 0
    while (max(abs(diff/y)) > eps & iter < mxiter) {
        sm = x + y + 1
        smd = sm^(-de1)
        G21 = (1 - smd)^(th1 - 1) * smd/sm
        gpdf = -G21
        gpdf = gpdf/(1 - smd)/sm/de/th
        gpdf = gpdf * (th * (de + 1) - (th * de + 1) * smd)
        iter = iter + 1
        diff = (G21 - pden)/gpdf
        y = y - diff
        while (min(y) <= 0 | max(abs(diff)) > 5) {
            diff = diff/2
            y = y + diff
        }
        if (iprint) 
            cat(iter, y, diff, "\n")
    }
    v = 1 - (1 - (y + 1)^(-de1))^th1
    if (iprint) 
        cat(v, "ended at iter. ", iter, "\n")
    v
}

qcondbb8 = function(p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) {
    vth = cpar[1]
    de = cpar[2]
    vth1 = 1/vth
    eta = (1 - (1 - de)^vth)
    eta1 = 1/eta
    be = 4 * pbb8(0.5, 0.5, cpar) - 1
    ut = (1 - de * u)^vth
    x = 1 - ut
    con = log(eta1) + (1 - vth1) * log(1 - x) - log(p)
    iter = 0
    diff = 1
    y = 0.5 * x
    if (be >= 0.8) 
        y = x
    while (iter < mxiter & abs(diff) > eps) {
        tem = 1 - eta1 * x * y
        h = log(y) + (vth1 - 1) * log(tem) + con
        hp = 1/y + (1 - vth1) * eta1 * x/tem
        diff = h/hp
        y = y - diff
        if (iprint) 
            cat(iter, diff, y, "\n")
        if (any(is.nan(y))) {
            if (iprint) 
                cat("***", p, u, cpar, x, "\n")
            return(u)
        }
        while (min(y) <= 0 | max(y) >= eta) {
            diff = diff/2
            y = y + diff
        }
        iter = iter + 1
    }
    if (iprint & iter >= mxiter) 
        cat("***did not converge\n")
    v = (1 - y)^vth1
    v = (1 - v)/de
    v
}

qcondbb9 = function(p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) {
    th = cpar[1]
    ga = cpar[2]
    ga1 = 1/ga
    x = ga1 - log(u)
    th1 = th - 1
    con = log(p) - x - th1 * log(x)
    z = (2 * x^th - ga1^th)^(1/th)
    mxdif = 1
    iter = 0
    diff = 0.1
    while (mxdif > eps & iter < mxiter) {
        h = z + th1 * log(z) + con
        hp = 1 + th1/z
        if (is.nan(h) | is.nan(hp) | is.nan(h/hp)) {
            diff = diff/-2
        }
        else diff = h/hp
        z = z - diff
        iter = iter + 1
        while (z <= x) {
            diff = diff/2
            z = z + diff
        }
        if (iprint) 
            cat(iter, diff, z, "\n")
        mxdif = max(abs(diff))
    }
    if (iprint) {
        if (iter >= mxiter) {
            cat("***did not converge ")
            cat("p, x, theta, ga, lastz :", p, x, th, ga, 
                z, "\n")
        }
    }
    y = (z^th - x^th + ga1^th)^(1/th)
    vv = exp(-y + ga1)
    if (iprint) 
        cat("p u v", p, u, vv, "\n")
    vv
}

qcondbb10 = function(p, u, cpar, eps = 1e-06, mxiter = 30, iprint = F) {
    th = cpar[1]
    ppi = cpar[2]
    th1 = 1/th
    be = 4 * pbb10(0.5, 0.5, cpar) - 1
    ut = u^th
    mxdif = 1
    iter = 0
    diff = 1
    v = 0.7 * u
    if (be >= 0.8) 
        v = u
    while (mxdif > eps & iter < mxiter) {
        vt = v^th
        tem = 1 - ppi * (1 - ut) * (1 - vt)
        ttem = tem^(-th1)
        ccdf = ttem/tem * v * (1 - ppi + ppi * vt)
        pdf = ttem/tem/tem
        pdf = pdf * (1 - ppi + ppi * (1 + th) * ut * vt - ppi * 
            (1 - ppi) * (1 - ut) * (1 - vt))
        h = ccdf - p
        hp = pdf
        diff = h/hp
        v = v - diff
        iter = iter + 1
        while (min(v) <= 0 | max(v) >= 1) {
            diff = diff/2
            v = v + diff
        }
        if (iprint) 
            cat(iter, diff, v, "\n")
        mxdif = abs(diff)
    }
    if (iprint & iter >= mxiter) {
        cat("***did not converge\n")
        cat("p=", p, " u=", u, " theta=", th, 
            " pi=", ppi, " lastv=", v, "\n")
    }
    v
}

qcondbvncop = function(p, v, cpar) {
    x = qnorm(p)
    y = qnorm(v)
    tem = x * sqrt(1 - cpar * cpar) + cpar * y
    pnorm(tem)
}

qcondbvtcop = function(p, u, cpar, df = dfdefault) {
    if (length(cpar) == 2) {
        rho = cpar[1]
        df = cpar[2]
    }
    else rho = cpar
    x2 = qt(p, df + 1)
    x1 = qt(u, df)
    tem = x2 * sqrt((df + x1 * x1) * (1 - rho * rho)/(df + 1)) + 
        rho * x1
    pt(tem, df)
}

qcondfrk = function(p, u, cpar) {
    cpar0 = exp(-cpar)
    cpar1 = 1 - cpar0
    etem = exp(-cpar * u + log(1/p - 1))
    tem = 1 - cpar1/(etem + 1)
    v = (-log(tem))/cpar
    isinf = is.infinite(v)
    v[isinf] = (-log(cpar0 + etem[isinf]))/cpar
    v
}

qcondjoe = function(p, u, cpar, eps = 1e-06, mxiter = 30) {
    cpar1 = cpar - 1
    cpartem = -cpar1/(1 + cpar1)
    cpar1inv = -1/cpar1
    ubar = 1 - u
    ud = ubar^cpar
    cpari = 1/cpar
    tem = (1 - p)^cpartem - 1
    tem = tem * (1 - u)^(-cpar1) + 1
    v = tem^cpar1inv
    v = 1 - v
    diff = 1
    iter = 0
    while (max(abs(diff)) > eps & iter < mxiter) {
        vbar = 1 - v
        vd = vbar^cpar
        sm = ud + vd - ud * vd
        smd = sm^cpari
        c21 = 1 + vd/ud - vd
        c21 = c21^(-1 + cpari)
        c21 = c21 * (1 - vd)
        pdf = smd * ud * vd * (cpar1 + sm)/sm/sm/ubar/vbar
        iter = iter + 1
        if (any(is.nan(pdf)) | any(is.nan(c21))) {
            diff = diff/(-2)
        }
        else diff = (c21 - p)/pdf
        v = v - diff
        while (min(v) <= 0 | max(v) >= 1 | max(abs(diff)) > 0.25) {
            diff = diff/2
            v = v + diff
        }
    }
    v
}

qcondmtcj = function(p, u, cpar) {
    eta = -cpar/(1 + cpar)
    tem = p^eta - 1
    tem = tem * (u^(-cpar)) + 1
    tem^(-1/cpar)
}

qcondpla = function(p, u, cpar, eps = 1e-08, mxiter = 30) {
    iter = 0
    diff = 1
    v = u
    while (iter < mxiter & max(abs(diff)) > eps) {
        num = pcondpla(v, u, cpar) - p
        den = dpla(u, v, cpar)
        diff = num/den
        v = v - diff
        while (min(v) < 0 | max(v) > 1) {
            diff = diff/2
            v = v + diff
        }
        iter = iter + 1
    }
    v
}

# bivariate copula conditional CDFs
pcondbb2 = function(v, u, cpar, iprint = F) {
    th = cpar[1]
    de = cpar[2]
    de1 = 1/de
    th1 = 1/th
    ut = u^(-th) - 1
    vt = v^(-th) - 1
    x = exp(ut * de) - 1
    y = exp(vt * de) - 1
    if (iprint) 
        cat(ut, vt, x, y, "\n")
    if (is.infinite(x) || is.infinite(y)) {
        lr = de * (vt - ut)
        r = exp(lr)
        tem = (1 + ut + de1 * log(1 + r))^(-th1 - 1)
        ccdf = tem * (ut + 1)/u/(1 + r)
    }
    else {
        sm = x + y + 1
        smd = de1 * log(sm)
        tem = (1 + smd)^(-th1 - 1)
        ccdf = tem * (x + 1) * (ut + 1)/sm/u
    }
    ccdf
}

pcondbb3 = function(v, u, cpar, iprint = F) {
    th = cpar[1]
    de = cpar[2]
    de1 = 1/de
    th1 = 1/th
    dt = (de^th1)
    ul = -log(u)
    vl = -log(v)
    ut = (ul^th)
    vt = (vl^th)
    x = exp(ut * de) - 1
    y = exp(vt * de) - 1
    if (iprint) 
        cat(ut, vt, x, y, "\n")
    if (is.infinite(x) || is.infinite(y) || x > 1e+200 || y > 
        1e+200) {
        lr = de * (vt - ut)
        r = exp(lr)
        lx = de * ut
        sm = r + 1
        sml = log(sm) + lx
        tem = (sml^th1)
        cdf = exp(-tem/dt)
        ccdf = cdf * tem * lx/sml/sm/ul/u/dt
    }
    else if (x <= 1e-10 || y <= 1e-10) {
        xx = de * ut
        yy = vt * de
        r = vt/ut
        sm = xx + yy
        tem = (1 + r)^(th1 - 1)
        cdf = exp(-(sm^th1)/dt)
        ccdf = cdf * tem * (1 + xx)/(1 + xx + yy)/u
    }
    else {
        sm = x + y + 1
        sml = log(sm)
        tem = (sml^th1)
        cdf = exp(-tem/dt)
        ccdf = cdf * tem * de * ut * (x + 1)/sml/sm/ul/u/dt
    }
    ccdf
}

pcondbb4 = function(v, u, cpar) {
    if (is.matrix(cpar)) {
        th = cpar[, 1]
        de = cpar[, 2]
    }
    else {
        th = cpar[1]
        de = cpar[2]
    }
    ut = u^(-th) - 1
    vt = v^(-th) - 1
    x = ut^(-de)
    y = vt^(-de)
    xy = x + y
    tem = xy^(-1/de)
    ccdf = ut + vt + 1 - tem
    xtem = ut/x - tem/xy
    ccdf = ccdf^(-1/th - 1) * xtem * x/ut * (1 + ut)/u
    ccdf
}

pcondbb5 = function(v, u, cpar) {
    if (is.matrix(cpar)) {
        th = cpar[, 1]
        de = cpar[, 2]
    }
    else {
        th = cpar[1]
        de = cpar[2]
    }
    x = -log(u)
    y = -log(v)
    xt = x^th
    yt = y^th
    xdt = x^(-de * th)
    ydt = y^(-de * th)
    xydt = xdt + ydt
    xyp = xydt^(-1/de)
    w = xt + yt - xyp
    wth = w^(1/th)
    cdf = exp(-wth)
    zx = xt/x - xyp/xydt * xdt/x
    ccdf = cdf * (wth/w) * zx/u
    ccdf
}

pcondbb6 = function(v, u, cpar) {
    if (is.matrix(cpar)) {
        th = cpar[, 1]
        de = cpar[, 2]
    }
    else {
        th = cpar[1]
        de = cpar[2]
    }
    ubar = 1 - u
    vbar = 1 - v
    zu = 1 - ubar^th
    x = -log(zu)
    y = -log(1 - vbar^th)
    xd = x^de
    yd = y^de
    sm = xd + yd
    tem = sm^(1/de)
    w = exp(-tem)
    ccdf = ((1 - w)/(1 - zu))^(1/th - 1) * (w/zu) * (tem/sm) * 
        (xd/x)
    ccdf
}

pcondbb7 = function(v, u, cpar) {
    if (is.matrix(cpar)) {
        th = cpar[, 1]
        de = cpar[, 2]
    }
    else {
        th = cpar[1]
        de = cpar[2]
    }
    de1 = 1/de
    th1 = 1/th
    ut = 1 - (1 - u)^th
    vt = 1 - (1 - v)^th
    x = ut^(-de) - 1
    y = vt^(-de) - 1
    sm = x + y + 1
    smd = sm^(-de1)
    tem = (1 - smd)^(th1 - 1)
    ccdf = tem * smd * (x + 1) * (1 - ut)/sm/ut/(1 - u)
    ccdf
}

pcondbb8 = function(v, u, cpar) {
    if (is.matrix(cpar)) {
        vth = cpar[, 1]
        de = cpar[, 2]
    }
    else {
        vth = cpar[1]
        de = cpar[2]
    }
    ut = (1 - de * u)^vth
    x = 1 - ut
    y = 1 - (1 - de * v)^vth
    eta1 = 1/(1 - (1 - de)^vth)
    tem = (1 - eta1 * x * y)^(1/vth - 1)
    den = (1 - de * u)/ut
    ccdf = eta1 * y * tem/den
    ccdf
}

pcondbb9 = function(v, u, cpar) {
    th = cpar[1]
    ga = cpar[2]
    ga1 = 1/ga
    x = ga1 - log(u)
    y = ga1 - log(v)
    temx = x^th
    temy = y^th
    sm = temx + temy - ga1^th
    smt = sm^(1/th)
    ccdf = exp(-smt + ga1)
    ccdf = ccdf * smt * temx/sm/x/u
    ccdf
}

pcondbb10 = function(v, u, cpar) {
    if (is.matrix(cpar)) {
        th = cpar[, 1]
        ppi = cpar[, 2]
    }
    else {
        th = cpar[1]
        ppi = cpar[2]
    }
    ut = u^th
    vt = v^th
    tem = 1 - ppi * (1 - ut) * (1 - vt)
    ttem = tem^(-1/th)
    ccdf = ttem/tem * v * (1 - ppi + ppi * vt)
    ccdf
}

pcondbvncop = function(v, u, cpar) {
    val = pnorm((qnorm(v) - cpar * qnorm(u))/sqrt(1 - cpar^2))
    val[v <= 0 | u <= 0 | u >= 1] = 0
    val[v == 1] = 1
    val
}

pcondbvtcop = function(v, u, cpar, df = dfdefault) {
    if (length(cpar) == 2) {
        rho = cpar[1]
        df = cpar[2]
    }
    else rho = cpar
    u[u == 0] = 1e-05
    v[v == 0] = 1e-06
    u[u == 1] = 1 - 1e-05
    v[v == 1] = 1 - 1e-06
    y1 = qt(u, df)
    y2 = qt(v, df)
    mu = rho * y1
    s2 = (1 - rho * rho) * (df + y1 * y1)/(df + 1)
    ccdf = pt((y2 - mu)/sqrt(s2), df + 1)
    ccdf
}

pcondfrk = function(v, u, cpar) {
    cpar[cpar == 0] = 1e-10
    cpar1 = 1 - exp(-cpar)
    tem = 1 - exp(-cpar * u)
    ccdf = (1 - tem)/(cpar1/(1 - exp(-cpar * v)) - tem)
    ccdf
}

pcondjoe = function(v, u, cpar) {
    temv = (1 - v)^cpar
    temu = (1 - u)^cpar
    ccdf = 1 + temv/temu - temv
    ccdf = ccdf^(-1 + 1/cpar)
    ccdf = ccdf * (1 - temv)
    ccdf
}

pcondmtcj = function(v, u, cpar) {
    tem = v^(-cpar) - 1
    tem = tem * (u^cpar) + 1
    ccdf = tem^(-1 - 1/cpar)
    ccdf
}

pcondpla = function (v, u, cpar) {
    cpar[cpar == 1] = 1 + 1e-10
    eta = cpar - 1
    tem = 1 + eta * (u + v)
    tem1 = tem * tem - 4 * cpar * eta * u * v
    tem2 = sqrt(tem1)
    ccdf = (eta * u + 1 - (eta + 2) * v)/tem2
    ccdf = 0.5 * (1 - ccdf)
    ccdf
}

################################################################################################################################################


# simulation for R-vine copulas when the first variable is set to a fixed value
# Arguments:
# p: vector of length d, e.g. runif(d)
# nsim: sample size for simulation
# A: d*d vine array with 1:d on diagonal, or ntrunc*d vine array as only ntrunc rows are used
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

rvinesimvec2_fix1 = function(nsim, A, ntrunc, parvec, np, qcondmat, pcondmat, extq, varname = numeric(0), iprint = F) {
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
