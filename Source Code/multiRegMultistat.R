# Function that fits multiple linear regression models simultaneously
# January 15, 2013

multiRegMultistat <- function (Ymat, Xmat, Llist, Penalty = NULL) 
{
    if (!is.matrix(Xmat)) {
        Xmat = cbind(1, Xmat)
        colnames(Xmat) = c("(Intercept)", "X")
    }
    nTest <- length(Llist)
    d = dim(Xmat)[2]
    if (is.null(colnames(Xmat))) 
        colnames(Xmat) = 1:d
    nobs = dim(Ymat)[2]
    nfeature = dim(Ymat)[1]
    Ymiss = is.na(Ymat)
    nmiss = apply(Ymiss, 1, sum)
    jallobs = which(nmiss == 0)
    janymiss = which(nmiss > 0)
    nanymiss = length(janymiss)
    beta = matrix(NA, nfeature, d)
    degfree = rep(NA, nfeature)
    stat = list()
    for (l in 1:nTest) stat[[l]] <- rep(NA, nfeature)
    colnames(beta) = colnames(Xmat)
    xtx = t(Xmat) %*% Xmat
    if (!is.null(Penalty)) 
        diag(xtx) = diag(xtx) + Penalty
    xtx.qr = qr(xtx)
    beta[jallobs, ] = t(solve(xtx.qr, t(Xmat) %*% t(Ymat[jallobs, 
        ])))
    Lbeta <- lapply(Llist, function(L) beta[jallobs, , drop = FALSE] %*% 
        L)
    LtXtXinvLqr <- lapply(Llist, function(L) qr(t(L) %*% solve(xtx.qr, 
        L)))
    for (l in 1:nTest) {
        for (j in 1:length(jallobs)) {
            stat[[l]][jallobs[j]] = Lbeta[[l]][j, ] %*% solve(LtXtXinvLqr[[l]], 
                t(Lbeta[[l]][j, , drop = FALSE]))
        }
    }
    degfree[jallobs] = nobs - d
    if (nanymiss > 0) 
        for (j in 1:nanymiss) {
            jj = janymiss[j]
            inc = which(!Ymiss[jj, ])
            if (length(inc) >= d) {
                XX2 = Xmat[inc, ]
                xtx = t(XX2) %*% XX2
                if (!is.null(Penalty)) 
                  diag(xtx) = diag(xtx) + Penalty
                xtx.qr = qr(xtx)
                beta[jj, ] = t(solve(xtx.qr, t(XX2) %*% as.numeric(Ymat[jj, 
                  inc])))
                Lbeta = sapply(Llist, function(L) beta[jj, , 
                  drop = FALSE] %*% L)
                for (l in 1:nTest) stat[[l]][jj] = Lbeta[[l]] %*% 
                  solve(t(Llist[[l]]) %*% solve(xtx.qr, Llist[[l]]), 
                    t(Lbeta[[l]]))
                degfree[jj] = length(inc) - d
            }
        }
    r = Ymat - beta %*% t(Xmat)
    sigmaSq = (1/degfree) * apply(r * r, 1, sum, na.rm = TRUE)
    standerrs = t(sapply(sigmaSq, function(w) sqrt(diag(solve(t(Xmat) %*% Xmat))*w)))
    stat = lapply(stat, function(st) st/sigmaSq)
    pvalue = sapply(stat, function(u) 1-(pf(u,1,degfree)))
    list(coef = beta, se = standerrs, stats = stat, degfree = degfree, pvalue = pvalue)
}