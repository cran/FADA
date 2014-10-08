decorrelate.train <- function (data.train, nbf = NULL, maxnbfactors = 12, nfolds = 10, 
    grouped = FALSE, plot.diagnostic = FALSE, min.err = 0.001, 
    verbose = TRUE,EM=TRUE,maxiter=15) 
{
	colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
                    twopass=FALSE) {
	  if (SumSquares) return(colSums(x^2, na.rm, dims))
	  N <- colSums(!is.na(x), FALSE, dims)
	  Nm1 <- if (unbiased) N-1 else N
	  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
	                     sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
	  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
	}
    bivprob = function(rho, lower, upper = -lower, mean = 0) {
        nu = 0
        low = rep(as.double((lower - mean)), 2)
        upp = rep(as.double((upper - mean)), 2)
        if (any(lower == upper)) 
            return(0)
        infin = c(2, 2)
        infin = as.integer(infin)
        low = replace(low, low == -Inf, 0)
        upp = replace(upp, upp == Inf, 0)
        rho = as.double(rho)
        prob = as.double(0)
        a = lapply(rho, function(r, low, upp) biv.nt.prob(df = Inf, 
            lower = low, upper = upp, mean = rep(0, 2), S = matrix(c(1, 
                r, r, 1), 2, 2)), low = low, upp = upp)
        return(unlist(a))
    }
    Dt = function(rho) {
        threshold = 0.05
        ut = qnorm(1 - threshold/2)
        delta = unlist(lapply(rho, bivprob, lower = -ut)) - (1 - 
            threshold)^2
        dt <- delta/(threshold * (1 - threshold))
        return(dt)
    }
    VarInflation <- function(data.train, Blist, maxnbfactors, dig) {
        m <- ncol(data.train)
        n <- nrow(data.train)
        vecrho <- round(seq(10^(-dig), 1, 10^(-dig)), digits = dig)
        vecdt <- unlist(lapply(vecrho, Dt))
        sampled <- sample(1:m, min(1000, m))
        sampsize <- length(sampled)
        cordata <- crossprod(data.train[, sampled, drop = FALSE])/(n - 
            1)
        sdt <- sapply(1:(maxnbfactors + 1), function(i) {
            B <- matrix(Blist[[i]][sampled, ], nrow = sampsize)
            sdb <- sqrt(1 - rowSums(B^2))  
            matrho <- cordata - tcrossprod(B)
            matrho <- sweep(matrho, 2, FUN = "/", STATS = sdb)
            matrho <- sweep(matrho, 1, FUN = "/", STATS = sdb)
            rho <- matrho[col(matrho) > row(matrho)]
            rho[abs(rho) >= 1] <- 1
            veccor <- sort(round(abs(rho), digits = dig))
            duplic <- duplicated(veccor)
            vduplic <- sort(unique(veccor[duplic]))
            vunic <- setdiff(unique(veccor), vduplic)
            dtunic <- vecdt[is.element(vecrho, vunic)]
            dtduplic <- vecdt[is.element(vecrho, vduplic)]
            vmatch <- match(vecrho, veccor, 0)
            nboccur <- diff(c(vmatch[vmatch > 0], length(veccor) + 
                1))
            nboccur <- nboccur[nboccur > 1]
            tmp <- 2 * (m - 1) * (sum(dtunic) + crossprod(nboccur, 
                dtduplic))/(sampsize * (sampsize - 1))
            return(tmp)
        })
        names(sdt) <- paste(0:maxnbfactors, "factors")
        return(sdt)
    }
    nbfactors <- function(data.train, maxnbfactors = 12, diagnostic.plot = plot.diagnostic, 
        minerr = 0.001, EM = TRUE) {
        dig <- 2
        m <- ncol(data.train)
        n <- nrow(data.train)
        falist <- vector(length = maxnbfactors + 1, "list")
        falist[[1]] <- list(B = matrix(0, ncol = 1, nrow = m))
        falist[-1] <- lapply(1:maxnbfactors, emfa, data = data.train, 
            minerr = minerr, EM = EM)
        Blist <- lapply(falist, function(fa, m) matrix(fa$B, 
            nrow = m), m = m)
        sdt <- VarInflation(data.train, Blist, maxnbfactors, dig)
        if (diagnostic.plot) {
            dev.new()
            plot(0:maxnbfactors, sdt, ylab = "Variance Inflation Criterion", 
                xlab = "Number of factors", bty = "l", lwd = 1.25, 
                type = "b", pch = 16, cex.lab = 1.25, cex = 1.25, 
                cex.axis = 1.25)
        }
        if (which.min(sdt) == 1) 
            opt <- 0
        if (which.min(sdt) > 1) {
            jumps <- -diff(sdt)/sdt[-length(sdt)]
            opt <- max((1:maxnbfactors)[jumps > 0.05])
        }
        list(criterion = sdt, optimalnbfactors = opt)
    }
    ifa = function(Psi, B) {
        if (class(B) == "numeric") 
            B = matrix(B, ncol = 1)
        q = ncol(B)
        Phi = rep(0, length(Psi))
        Phi[abs(Psi) > 1e-05] = 1/Psi[abs(Psi) > 1e-05]
        PhiB = tcrossprod(Phi,rep(1, q))
        PhiB = PhiB * B
        G = diag(q) + t(B) %*% PhiB
        GinvtPhiB = tcrossprod(solve(G),PhiB)
        Phib2 = tcrossprod(PhiB, t(GinvtPhiB))
        iS = diag(Phi) - Phib2
        PhiB2 = crossprod(PhiB, B)
        GinvtPhiB2 = crossprod(solve(G),PhiB2) 
        Phib2 = tcrossprod(PhiB, t(GinvtPhiB2))
        iSB = PhiB - Phib2
        return(list(iS = iS, iSB = iSB))
    }
    emfa = function(data, nbf, EM = TRUE, minerr = 1e-06, verbose = FALSE) {
        n = nrow(data)
        m = ncol(data)
        my = crossprod(rep(1, n), data)/n
        vy = crossprod(rep(1, n), data^2)/n - my^2
        vy = (n/(n - 1)) * vy
        cdata = scale(data, center = my, scale = FALSE)
        csdata = scale(data, center = my, scale = sqrt(vy))
        S = crossprod(csdata)/(n - 1)
        if (((n > m) & (m <= 200) & (m >= 3)) & (!EM)) {
            if (nbf == 0) {
                B = NULL
                Psi = rep(1, m)
            }
            if (nbf > 0) {
                fa = factanal(csdata, factors = nbf, rotation = "varimax")
                B = fa$loadings
                class(B) = "matrix"
                Psi = fa$uniquenesses
                Psi = Psi*vy
                B = matrix(rep(sqrt(vy),ncol(B)),nrow=nrow(B))*B
                sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
                G = solve(diag(nbf) + tcrossprod(sB))
                sB = scale(t(B), center = FALSE, scale = Psi)
            }
        }
        if ((n <= m) | (m > 200) | EM) {
            if (nbf == 0) {
                B = NULL
                Psi = rep(1, m)
            }
            if (nbf > 0) {
                if (verbose) 
                  print(paste("Fitting EM Factor Analysis Model with", 
                    nbf, "factors"))
                eig = svd((1/sqrt((n - 1))) * t(csdata))
                evectors = eig$u[, 1:nbf]
                evalues = eig$d^2
                if (nbf > 1) 
                  B = evectors[, 1:nbf] %*% diag(sqrt(evalues[1:nbf]))
                if (nbf == 1) 
                  B = matrix(evectors, nrow = m, ncol = 1) * 
                    sqrt(evalues[1])
                b2 = rowSums(B^2) 
                Psi = 1 - b2
                crit = 1
                while (crit > minerr) {
                  inv = ifa(Psi, B)
                  Cyz = crossprod(S ,inv$iSB)
                  Czz = crossprod(inv$iSB, Cyz) + diag(nbf) - 
                    crossprod(B, inv$iSB)
                  Bnew = tcrossprod(Cyz ,solve(Czz))
                  Psinew = 1 - rowSums(Bnew * Cyz) 
                  crit = mean((Psi - Psinew)^2)
                  B = Bnew
                  Psi = Psinew
                  if (verbose) 
                    print(paste("Objective criterion in EM-FA : ", 
                      signif(crit, 6)))
                }
                Psi = Psi*vy
                B = matrix(rep(sqrt(vy),ncol(B)),nrow=nrow(B))*B
                sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
                G = solve(diag(nbf) + tcrossprod(sB))
                sB = scale(t(B), center = FALSE, scale = Psi)
            }
        }
       res = list(B = B, Psi = Psi)
        return(res)
    }
    LassoML <- function(data.train, type.multinomial = "ungrouped",nfold=nfolds) {
        p <- ncol(data.train$x)
        n <- nrow(data.train$x)
        nbclass <- length(unique(data.train$y))
        cl <- sort(unique(data.train$y))
        if (!all(cl == c(1:nbclass))) {
            stop
            ("Group variable must be 1,2, ...")
        }
        family <- ifelse(nbclass==2,"binomial","multinomial")
            cvmod <- cv.glmnet(x = as.matrix(data.train$x), y = data.train$y, 
                family = family, type.measure = "class", 
                nfolds = nfolds, grouped = FALSE,type.multinomial = type.multinomial)
            lambda.min <- cvmod$lambda.min
            mod <- glmnet(x = as.matrix(data.train$x), y = data.train$y, family = family, 
                lambda = lambda.min,type.multinomial = type.multinomial)
            proba.train <- predict(mod, newx = as.matrix(data.train$x), 
                type = "response")
       if (nbclass == 2) { proba.train <- matrix(c(1 - proba.train, proba.train), ncol = 2, byrow = FALSE)} 
       if (nbclass > 2) {proba.train <- proba.train[, , 1]}
        return(list(proba.train = proba.train,model = mod))
    }
    p <- ncol(data.train$x)
    n <- nrow(data.train$x)
    nbclass <- length(unique(data.train$y))
    cl <- sort(unique(data.train$y))
    if (!all(cl == c(1:nbclass))) {
        stop("Group variable must be 1,2, ...")
    }
    if (p <= 3) {stop("The number of variables must be at least 4.")}
    meanclass <- sapply(1:nbclass, function(i) {
        colMeans(data.train$x[data.train$y == i, ])
    })
    cdta <- data.train$x
    for (i in 1:nbclass) {
        cdta[data.train$y == i, ] <- sweep(data.train$x[data.train$y == i, ], 2, meanclass[, 
            i], "-")
    }
    if (is.null(nbf)) {
        nbf <- nbfactors(scale(cdta, center = FALSE, scale = TRUE), 
            maxnbfactors = maxnbfactors, EM = EM,diagnostic.plot=plot.diagnostic)$optimalnbfactors
    }
    if (verbose) 
        print(paste("Number of factors:", nbf, "factors", sep = " "))
    coefpx <- LassoML(data.train,nfold=nfolds)
    mod <- coefpx$mod
    if (nbf == 0) {    	
        return(list(meanclass = meanclass, fa.training = data.train$x, 
            Psi = colVars(cdta),
            B = NULL,factors.training=NULL, groups = data.train$y,proba.training=coefpx$proba.train,mod.decorrelate.test=mod,data.train=data.train$x))
    }
    else {
        eps <- 1
        iter <- 0
        proba.train <- coefpx$proba.train
        if (verbose) 
            print("Objective criterion: ")
        while (eps > min.err & iter <= maxiter) {
            fa <- emfa(cdta, nbf = nbf, minerr = min.err,verbose=FALSE,EM=EM)
            sB <- scale(t(fa$B), center = FALSE, scale = sqrt(fa$Psi))
            G <- solve(diag(ncol(fa$B)) + tcrossprod(sB))
            sB <- scale(t(fa$B), center = FALSE, scale = fa$Psi)
            zclass <- lapply(1:nbclass, function(i) {
                sweep(tcrossprod(tcrossprod(scale(data.train$x, center = meanclass[, 
                  i], scale = FALSE), sB), G), 1, proba.train[, 
                  i], "*")
            })
            z <- Reduce("+", zclass)
            fadta <- data.train
            fadta$x <- data.train$x - tcrossprod(z, fa$B)
            fameanclass <- sapply(1:nbclass, function(i) {
                colMeans(fadta$x[data.train$y == i, ])
            })
            coefpx <- LassoML(fadta,nfold=nfolds)
            proba.train <- coefpx$proba.train
            v = colVars(cdta)
            eps <- colSums((fameanclass/sqrt(v)-meanclass/sqrt(v))^2)
            eps <-(mean(eps))
            iter <- iter+1
            if (verbose) 
                print(eps)
            meanclass <- fameanclass
            cdta <- data.train$x
            for (i in 1:nbclass) {
                cdta[data.train$y == i, ] <- sweep(data.train$x[data.train$y == i, 
                  ], 2, meanclass[, i], "-")
            } 
        }
        fa <- emfa(cdta, nbf = nbf, minerr = 1e-06,EM=EM)
        sB <- scale(t(fa$B), center = FALSE, scale = sqrt(fa$Psi))
        G <- solve(diag(ncol(fa$B)) + tcrossprod(sB))
        sB <- scale(t(fa$B), center = FALSE, scale = fa$Psi)
        zclass <- lapply(1:nbclass, function(i) {
            tcrossprod(tcrossprod(scale(data.train$x, center = meanclass[, 
                i], scale = FALSE), sB), G)
        })
        zclass <- lapply(1:nbclass, function(i) {
            sweep(zclass[[i]], 1, proba.train[, i], "*")
        })
        z <- Reduce("+", zclass)
        fadta <- data.train
        fadta$x <- data.train$x - tcrossprod(z, fa$B)
        return(list(meanclass = meanclass, fa.training = fadta$x ,
            Psi = fa$Psi, B = fa$B,factors.training=z,groups = data.train$y,proba.training= proba.train,mod.decorrelate.test=mod,data.train=data.train$x))
    }
}
