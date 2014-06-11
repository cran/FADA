FA = function(dta,test=NULL,nbf=NULL,maxnbfactors=12,nfolds=10,grouped=FALSE,plot.diagnostic=FALSE,min.err = 0.001,verbose=TRUE){
	bivprob = function (rho,lower,upper=-lower,mean=0) {
    nu = 0
    low = rep(as.double((lower - mean)),2)
    upp = rep(as.double((upper - mean)),2)
    if (any(lower == upper)) return(0)
    infin = c(2, 2)
    infin = as.integer(infin)
    low = replace(low, low == -Inf, 0)
    upp = replace(upp, upp == Inf, 0)
    rho = as.double(rho)
    prob = as.double(0)
    a = lapply(rho,function(r,low,upp) biv.nt.prob(df=Inf,lower=low,upper=upp,mean=rep(0,2),S=matrix(c(1,r,r,1),2,2)),
               low=low,upp=upp)
    return(unlist(a))
  }
  
  Dt = function(rho) {
    threshold=0.05
    ut = qnorm(1 - threshold/2)
    delta = unlist(lapply(rho,bivprob,lower=-ut)) - (1 - threshold)^2
    dt <-  delta/(threshold * (1 - threshold))
    return(dt)
  }
  
  VarInflation <-  function(dta,Blist,maxnbfactors,dig) {
    m <-  ncol(dta)
    n <-  nrow(dta)
    vecrho <-  round(seq(10^(-dig),1,10^(-dig)),digits=dig)
    vecdt <-  unlist(lapply(vecrho,Dt))
    sampled <-  sample(1:m,min(1000,m))
    sampsize <-  length(sampled)
    cordata <-  crossprod(dta[,sampled,drop=FALSE])/(n-1)
    sdt <-  sapply(1:(maxnbfactors+1),function(i) {
      B <-  matrix(Blist[[i]][sampled,],nrow=sampsize)
      sdb <-  sqrt(1-apply(B^2,1,sum))
      matrho <-  cordata - tcrossprod(B)
      matrho <-  sweep(matrho,2,FUN="/",STATS=sdb)
      matrho <-  sweep(matrho,1,FUN="/",STATS=sdb)
      rho <-  matrho[col(matrho)>row(matrho)]
      rho[abs(rho)>=1] <-  1
      veccor <-  sort(round(abs(rho),digits=dig))
      duplic <-  duplicated(veccor)
      vduplic <-  sort(unique(veccor[duplic]))
      vunic <-  setdiff(unique(veccor),vduplic)
      dtunic <-  vecdt[is.element(vecrho,vunic)]
      dtduplic <-  vecdt[is.element(vecrho,vduplic)]
      vmatch <-  match(vecrho,veccor,0)
      nboccur <-  diff(c(vmatch[vmatch>0],length(veccor)+1))
      nboccur <-  nboccur[nboccur>1]
      tmp <-  2*(m-1)*(sum(dtunic)+crossprod(nboccur,dtduplic))/(sampsize*(sampsize-1))  
      return(tmp)})
      names(sdt) <-  paste(0:maxnbfactors,"factors")
      return(sdt)
     }
  
  nbfactors <-  function(dta,maxnbfactors=12,diagnostic.plot=plot.diagnostic,minerr=1e-03){
    dig <-  2
    m <-  ncol(dta)
    n <-  nrow(dta)
    falist <-  vector(length=maxnbfactors+1,"list")
    falist[[1]] <-  list(B=matrix(0,ncol=1,nrow=m))
    falist[-1] <-  lapply(1:maxnbfactors,emfa,data=dta,minerr=minerr)
    Blist <-  lapply(falist,function(fa,m) matrix(fa$B,nrow=m),m=m)
    sdt <-  VarInflation(dta,Blist,maxnbfactors,dig)
    if (diagnostic.plot) {
      dev.new()
      plot(0:maxnbfactors,sdt,ylab="Variance Inflation Criterion",xlab="Number of factors",bty="l",
           lwd=1.25,type="b",pch=16,cex.lab=1.25,cex=1.25,cex.axis=1.25)
    }
    if (which.min(sdt)==1) opt <-  0
    if (which.min(sdt)>1) { 
      jumps <-  -diff(sdt)/sdt[-length(sdt)]
      opt <-  max((1:maxnbfactors)[jumps > 0.05])
    }
    list(criterion=sdt,optimalnbfactors=opt)
  }
  

ifa = function(Psi,B) {
  if (class(B)=="numeric") B = matrix(B,ncol=1)
  q = ncol(B)
  Phi = rep(0,length(Psi))
  Phi[abs(Psi)>1e-05] = 1/Psi[abs(Psi)>1e-05]
  PhiB = Phi%*%t(rep(1,q))
  PhiB = PhiB*B
  G = diag(q)+t(B)%*%PhiB
  GinvtPhiB = solve(G)%*%t(PhiB)
  Phib2 = tcrossprod(PhiB,t(GinvtPhiB))
  iS = diag(Phi) - Phib2
  PhiB2 = crossprod(PhiB,B)
  GinvtPhiB2 = solve(G)%*%PhiB2
  Phib2 = tcrossprod(PhiB,t(GinvtPhiB2))
  iSB = PhiB - Phib2
  return(list(iS=iS,iSB=iSB))
}

emfa = function (data, nbf, EM=TRUE,minerr = 1e-06,verbose=FALSE) {
  n = nrow(data)
  m = ncol(data)
  my = crossprod(rep(1,n),data)/n
  vy = crossprod(rep(1,n),data^2)/n-my^2
  vy = (n/(n-1))*vy
  cdata = scale(data,center=my,scale=FALSE)
  S = crossprod(cdata)/(n-1)
  if ((n>m)&(m<=200)&(m>=3)&(!EM)) {
    if (nbf == 0) {
      B = NULL
      Psi = rep(1, m)
      Factors = NULL
    }
    if (nbf > 0) {
      #           print(paste("Fitting direct ML Factor Analysis Model with", nbf,"factors"))
      fa = factanal(data,factors=nbf,rotation="varimax")
      B = fa$loadings
      class(B) = "matrix"
      Psi = fa$uniquenesses
      sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
      G = solve(diag(nbf) + sB %*% t(sB))
      sB = scale(t(B), center = FALSE, scale = Psi)
      Factors = cdata%*%t(sB)%*%t(G)
    }   
  }
  if ((n<=m)|(m>200)|(m<=2)|EM) {
    if (nbf == 0) {
      B = NULL
      Psi = rep(1, m)
      Factors = NULL
    }
    if (nbf > 0) {
      if (verbose) print(paste("Fitting EM Factor Analysis Model with", nbf,"factors"))
      eig = svd((1/sqrt((n - 1))) * t(cdata))
      evectors = eig$u[, 1:nbf]
      evalues = eig$d^2
      if (nbf > 1) 
        B = evectors[, 1:nbf] %*% diag(sqrt(evalues[1:nbf]))
      if (nbf == 1) 
        B = matrix(evectors, nrow = m, ncol = 1) * sqrt(evalues[1])
      b2 = apply(B^2, 1, sum)
      Psi = vy - b2
      crit = 1
      while (crit > minerr) {
        inv = ifa(Psi,B)
        Cyz = S%*%inv$iSB
        Czz = crossprod(inv$iSB,Cyz)+diag(nbf)-crossprod(B,inv$iSB)
        Bnew = Cyz%*%solve(Czz)
        Psinew = vy - apply(Bnew*Cyz,1,sum)
        crit = mean((Psi - Psinew)^2)
        B = Bnew
        Psi = Psinew
        if (verbose) print(paste("Objective criterion in EM-FA : ",signif(crit,6)))
      }
      sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
      G = solve(diag(nbf) + sB %*% t(sB))
      sB = scale(t(B), center = FALSE, scale = Psi)
      Factors = cdata%*%t(sB)%*%t(G)
    }
  }
  res = list(B = B, Psi = Psi,Factors=Factors)
  return(res)
}

    LassoML <-  function(dta,test=NULL,type.multinomial="ungrouped") {
  p <-  ncol(dta$x)
  n <-  nrow(dta$x)
  nbclass <-  length(unique(dta$y))
  cl <-  sort(unique(dta$y))
  if (! all(cl == c(1:nbclass))) {stop
                                  ("Group variable must be 1,2, ...")}
  if (nbclass == 2) {
	    cvmod <-  cv.glmnet(x=as.matrix(dta$x),y=dta$y,family="binomial",type.measure="class",nfolds=n,grouped=FALSE)
	    lambda.min <-  cvmod$lambda.min
	    mod <-  glmnet(x=as.matrix(dta$x),y=dta$y,family="binomial",lambda=lambda.min)
	    proba.train <-  predict(mod,newx=as.matrix(dta$x),type="response")
	    proba.train <-  matrix(c(1-proba.train,proba.train),ncol=2,byrow=FALSE)
	    if (is.null(test)) {proba.test=NULL;predict.test=NULL} 
	    else {
	    proba.test <-  predict(mod,newx=as.matrix(test$x),type="response")
	    proba.test <-  matrix(c(1-proba.test, proba.test),ncol=2,byrow=FALSE)
	    predict.test <-  apply(proba.test,1,which.max)
	    }
	    selected <-  which(abs(coef(mod)[-1])>1e-06) 
	    }
	  else { 
	  	cvmod <-  cv.glmnet(x=as.matrix(dta$x),y=dta$y,family="multinomial",type.measure="class",type.multinomial=type.multinomial,nfolds=n,grouped=FALSE)
	    lambda.min <-  cvmod$lambda.min
	    mod <-  glmnet(x=as.matrix(dta$x),y=dta$y,family="multinomial",lambda=lambda.min,type.multinomial=type.multinomial)
	    proba.train <-  predict(mod,newx=as.matrix(dta$x),type="response")[,,1]
	    if (is.null(test)) {proba.test<- NULL;predict.test<- NULL} 
	    else {
	    proba.test <-  predict(mod,newx=as.matrix(test$x),type="response")[,,1]
	    if (nrow(test$x) == 1) {proba.test <-  matrix(proba.test,nrow=1)}
	    predict.test <-  apply(proba.test,1,which.max)
	    }
	    selected <-  lapply(coef(mod),function(x) which(abs(x[-1])>1e-06 ) ) 
	    } ; 	    return(list(proba.train=proba.train,proba.test=proba.test,predict.test=predict.test,selected=selected,model=mod))
}
  
  p <-  ncol(dta$x)
  n <-  nrow(dta$x)
  nbclass <-  length(unique(dta$y))
  cl <-  sort(unique(dta$y))
  if (! all(cl == c(1:nbclass))) {stop("Group variable must be 1,2, ...")}
  
  meanclass <-  sapply(1:nbclass,function(i){apply(dta$x[dta$y==i,],2,mean)})
  
  cdta <-  dta$x
  for (i in 1:nbclass){
    cdta[dta$y==i,] <-  sweep(dta$x[dta$y==i,],2,meanclass[,i],"-")  
  }
  
  if (is.null(nbf)) { nbf <-  nbfactors(scale(cdta,center=FALSE,scale=TRUE),maxnbfactors=maxnbfactors)$optimalnbfactors }
  if (verbose) print(paste("Number of factors:",nbf,"factors",sep=" "))
  
  coefpx <-  LassoML(dta=dta,test=test) 
  proba.test <-  coefpx$proba.test
  predict.test <-  coefpx$predict.test
  
  if (nbf == 0) { if (is.null(test)) {dta.tmp<- dta$x} else {dta.tmp <-  NULL}
    return(list(meanclass=meanclass,
                proba.test=proba.test,predict.test=predict.test,fadta=dta$x,fatest=test$x,Psi=apply(cdta,2,var),B=NULL,groups=dta$y,dta=dta.tmp))
  }
  else {
    eps <-  1
    proba.train <-  coefpx$proba.train
    
    if (verbose) print("Objective criterion: ")
    while (eps>min.err) {
      fa <-  emfa(cdta,nbf=nbf,minerr=1e-06)
      
      sB <-  scale(t(fa$B), center = FALSE, scale = sqrt(fa$Psi))
      G  <-  solve(diag(ncol(fa$B)) + tcrossprod(sB))
      sB <-  scale(t(fa$B), center = FALSE, scale = fa$Psi)
      
      zclass <-  lapply(1:nbclass,function(i){sweep(tcrossprod(tcrossprod(scale(dta$x,center=meanclass[,i],scale=FALSE),sB),G),1,proba.train[,i],"*")})
      z <-  Reduce('+', zclass)
      
      fadta <-  dta
      fadta$x <-  dta$x-tcrossprod(z,fa$B)       
      
      fameanclass <-  sapply(1:nbclass,function(i){apply(fadta$x[dta$y==i,],2,mean)})
      
      coefpx <-  LassoML(dta=fadta,test=test) 
      
      proba.train <-  predict(coefpx$mod,dta$x,type="response")    	
      if (nbclass == 2) {proba.train <-  matrix(c(1-proba.train,proba.train),ncol=2,byrow=FALSE)}
      proba.train <-  matrix(proba.train,ncol=nbclass,byrow=FALSE) 
      
      eps <-  apply((fameanclass - meanclass)^2,1,sum)
      eps <-  sqrt(mean(eps))
      if (verbose) print(eps)
      
      meanclass <-  fameanclass
      
      cdta <-  dta$x
      for (i in 1:nbclass){
        cdta[dta$y==i,] <-  sweep(dta$x[dta$y==i,],2,meanclass[,i],"-")  
      }    
    }
    
    fa <-  emfa(cdta,nbf=nbf,minerr=1e-06)
    
    sB <-  scale(t(fa$B), center = FALSE, scale = sqrt(fa$Psi))
    G  <-  solve(diag(ncol(fa$B)) + tcrossprod(sB))
    sB <-  scale(t(fa$B), center = FALSE, scale = fa$Psi)
    
    zclass <-  lapply(1:nbclass,function(i){tcrossprod(tcrossprod(scale(dta$x,center=meanclass[,i],scale=FALSE),sB),G)})
    zclass <-  lapply(1:nbclass,function(i){sweep(zclass[[i]],1,proba.train[,i],"*")})
    z <-  Reduce('+', zclass)
    
    fadta <-  dta
    fadta$x <-  dta$x-tcrossprod(z,fa$B)       
    
    if (is.null(test)) {fatest<- NULL} else {
    testzclass <-  lapply(1:nbclass,function(i){tcrossprod(tcrossprod(scale(test$x,center=meanclass[,i],scale=FALSE),sB),G)})
    testzclass <-  lapply(1:nbclass,function(i){sweep(testzclass[[i]],1,proba.test[,i],"*")})
    testz <-  Reduce('+', testzclass)
    
    fatest <-  test
    fatest$x <-  test$x-tcrossprod(testz,fa$B) 
    
    }
    if (is.null(test)) {dta.tmp<- dta$x} else {dta.tmp<- NULL}
    return(list(meanclass=meanclass,fadta=fadta$x,fatest=fatest$x,Psi=fa$Psi,B=fa$B,groups=dta$y,dta=dta.tmp))	
  }
}



FADA <-  function(faobject,nfold.cv, nbf.cv=NULL, method = c("glmnet", "sda", "sparseLDA"), 
    sda.method = c("lfdr", "HC"), stop.par = 10, lambda, lambda.var, 
    lambda.freqs, diagonal = FALSE, alpha = 0.1,nfolds=10){
  
  nbclass <- length(unique(faobject$groups))
      	
         LassoML <-  function(dta,test=NULL,type.multinomial="ungrouped") {
  p <-  ncol(dta$x)
  n <-  nrow(dta$x)
  nbclass <-  length(unique(dta$y))
  cl <-  sort(unique(dta$y))
  if (! all(cl == c(1:nbclass))) {stop
                                  ("Group variable must be 1,2, ...")}
  if (nbclass == 2) {
	    cvmod <-  cv.glmnet(x=as.matrix(dta$x),y=dta$y,family="binomial",type.measure="class",nfolds=n,grouped=FALSE)
	    lambda.min <-  cvmod$lambda.min
	    mod <-  glmnet(x=as.matrix(dta$x),y=dta$y,family="binomial",lambda=lambda.min)
	    proba.train <-  predict(mod,newx=as.matrix(dta$x),type="response")
	    proba.train <-  matrix(c(1-proba.train,proba.train),ncol=2,byrow=FALSE)
	    if (is.null(test)) {proba.test=NULL;predict.test=NULL} 
	    else {
	    proba.test <-  predict(mod,newx=as.matrix(test$x),type="response")
	    proba.test <-  matrix(c(1-proba.test, proba.test),ncol=2,byrow=FALSE)
	    predict.test <-  apply(proba.test,1,which.max)
	    }
	    selected <-  which(abs(coef(mod)[-1])>1e-06) 
	    }
	  else { 
	  	cvmod <-  cv.glmnet(x=as.matrix(dta$x),y=dta$y,family="multinomial",type.measure="class",type.multinomial=type.multinomial,nfolds=n,grouped=FALSE)
	    lambda.min <-  cvmod$lambda.min
	    mod <-  glmnet(x=as.matrix(dta$x),y=dta$y,family="multinomial",lambda=lambda.min,type.multinomial=type.multinomial)
	    proba.train <-  predict(mod,newx=as.matrix(dta$x),type="response")[,,1]
	    if (is.null(test)) {proba.test<- NULL;predict.test<- NULL} 
	    else {
	    proba.test <-  predict(mod,newx=as.matrix(test$x),type="response")[,,1]
	    if (nrow(test$x) == 1) {proba.test <-  matrix(proba.test,nrow=1)}
	    predict.test <-  apply(proba.test,1,which.max)
	    }
	    selected <-  lapply(coef(mod),function(x) which(abs(x[-1])>1e-06 ) ) 
	    } ; 	    return(list(proba.train=proba.train,proba.test=proba.test,predict.test=predict.test,selected=selected,model=mod))
}

    	
FADA.tmp <-  function (faobject, method, sda.method, stop.par , lambda, lambda.var,  lambda.freqs, diagonal , alpha ) 
{
    fadta <-  faobject$fadta
    fatest <-  faobject$fatest
    groups <-  faobject$groups
    p <-  ncol(faobject$fadta)
    if (method == "glmnet") {
        out <-  LassoML(list(x = fadta, y = groups), list(x = fatest))
        selected <-  out$selected
        proba.test <-  out$proba.test
        predict.test <-  out$predict.test
        beta0 <-  coef(out$model)[1]
        beta <-  coef(out$model)[-1]
        beta <- lapply(1:(nbclass-1),function(k){beta[[k]][selected[[k]]]})
    }
    if (method == "sda") {
        ranking.LDA <-  sda::sda.ranking(fadta, groups, diagonal = diagonal, 
            verbose = FALSE)
        if (sda.method == "lfdr") {
            selected <-  as.numeric(ranking.LDA[ranking.LDA[, "lfdr"] < 
                0.8, "idx"])
        }
        else {
            thr <-  which.max(ranking.LDA[1:round(alpha * p), "HC"])
            selected <-  as.numeric(ranking.LDA[1:thr, "idx"])
        }
        out <-  sda::sda(fadta[, selected], groups, verbose = FALSE)
        pred <-  sda::predict.sda(out, fatest[, selected,drop=FALSE], verbose = FALSE)
        proba.test <-  pred$posterior
        predict.test <-  pred$class
        beta0 <-  out$alpha
        beta <-  out$beta
    }
    if (method == "sparseLDA") {
        Xc <-  normalize(fadta)
        Xn <- Xc$Xc
        out <-  sda(Xn, factor(groups), stop = -stop.par)
        Xctest <-  normalizetest(fatest, Xc)
        Xctest <-  matrix(Xctest, nrow = nrow(fatest), byrow = FALSE)
        colnames(Xctest) <-  colnames(Xn)
        pred <-  predict.sda(out, Xctest)
        selected <-  out$varIndex
        proba.test <-  pred$posterior
        predict.test <-  pred$class
        beta0 <-  NULL
        beta <-  out$beta
    }
    return(list(method = method, selected = selected, proba.test = proba.test, 
        predict.test = predict.test, beta0 = beta0, beta = beta))
}
    	

cv.FADA <-  function (dta, nfold.cv= nfold.cv, nbf.cv=nbf.cv, method = c("glmnet", 
    "sda", "sparseLDA"), sda.method = c("lfdr", "HC"), stop.par = 10, 
    nfolds = 10, lambda, lambda.var, lambda.freqs, 
    diagonal = FALSE, alpha = 0.1) 
{
    n <-  nrow(dta$x)
    p <-  ncol(dta$x)
    if (nfold.cv == n) {
        fold <-  as.list(c(1:n))
    }
    else {
        fold <-  lapply(c(1:nfold.cv), function(j) {
            sample(x = c(1:n), size = ceiling(n/2))
        })
        fold <-  lapply(fold, sort)
    }
    dta.cv <-  list()
    test.cv <-  list()
    cv.func <- function(j) {
    	print(j)
        dta.cv$x <-  dta$x[-fold[[j]], ]
        dta.cv$y <-  dta$y[-fold[[j]]]
        test.cv$x <-  dta$x[fold[[j]], ,drop=FALSE]
        fa.cv <-  FA(dta.cv, test.cv, nbf = nbf.cv, nfolds = nfolds, 
            verbose = FALSE)
        test.cv$y <-  dta$y[fold[[j]]]
        fada <-  FADA.tmp(fa.cv, method = method, sda.method = sda.method, 
            stop.par = stop.par, lambda = lambda, lambda.var = lambda.var, 
            lambda.freqs = lambda.freqs, diagonal = diagonal, 
            alpha = alpha)
        return(mean(fada$predict.test != test.cv$y))
    }
    cv.fa <- sapply(1:nfold.cv,cv.func)
    return(mean(cv.fa))
}

    if (is.null(faobject$dta)) {out <-  FADA.tmp(faobject, method,  sda.method , stop.par , lambda, lambda.var, lambda.freqs, diagonal, alpha)
    		proba.test <-  out$proba.test
			predict.test <-  out$predict.test
			selected <-  out$selected
			beta0 <-  out$beta0
			beta <-  out$beta
    		cv.error <-  NULL }
    	else { 
    		fadta <-  faobject$fadta
    		groups <-  faobject$groups
		    p <-  ncol(faobject$fadta)
		    if (method == "glmnet") {
		        out <-  LassoML(list(x = fadta, y = groups), NULL)
		        selected <-  out$selected
		        beta0 <-  coef(out$model)[1]
		        beta <-  coef(out$model)[-1]
		        beta <- lapply(1:(nbclass-1),function(k){beta[[k]][selected[[k]]]})
		    }
		    if (method == "sda") {
		        ranking.LDA <-  sda::sda.ranking(fadta, groups, diagonal = diagonal, 
		            verbose = FALSE)
		        if (sda.method == "lfdr") {
		            selected <-  as.numeric(ranking.LDA[ranking.LDA[, "lfdr"] < 
		                0.8, "idx"])
		        }
		        else {
		            thr <-  which.max(ranking.LDA[1:round(alpha * p), "HC"])
		            selected <-  as.numeric(ranking.LDA[1:thr, "idx"])
		        }
		        out <-  sda::sda(fadta[, selected], groups, verbose = FALSE)
		        beta0 <-  out$alpha
		        beta <-  out$beta
		    }
		    if (method == "sparseLDA") {
		        Xc <-  normalize(fadta)
		        Xn <- Xc$Xc
		        out <-  sda(Xn, factor(groups), stop = -stop.par)
		        selected <-  out$varIndex
		        beta0 <-  NULL
		        beta <-  out$beta
		    }
	cv.out <-  cv.FADA(dta=list(x=faobject$dta,y=faobject$groups),nfold.cv=nfold.cv, nbf.cv=nbf.cv, method=method, sda.method=sda.method, stop.par= stop.par ,  nfolds= nfolds, lambda= lambda, lambda.var= lambda.var, lambda.freqs= lambda.freqs, diagonal=diagonal, alpha=alpha)
	proba.test <-  NULL
	predict.test <-  NULL
	cv.error <-  cv.out
    	}
return(list(method = method, selected = selected, proba.test = proba.test, 
        predict.test = predict.test,cv.error=cv.error, beta0 = beta0, beta = beta))
    }
    
    
decorrelate = function(faobject,test){
	   LassoML <-  function(dta,test=NULL,type.multinomial="ungrouped") {
  p <-  ncol(dta$x)
  n <-  nrow(dta$x)
  nbclass <-  length(unique(dta$y))
  cl <-  sort(unique(dta$y))
  if (! all(cl == c(1:nbclass))) {stop
                                  ("Group variable must be 1,2, ...")}
  if (nbclass == 2) {
	    cvmod <-  cv.glmnet(x=as.matrix(dta$x),y=dta$y,family="binomial",type.measure="class",nfolds=n,grouped=FALSE)
	    lambda.min <-  cvmod$lambda.min
	    mod <-  glmnet(x=as.matrix(dta$x),y=dta$y,family="binomial",lambda=lambda.min)
	    proba.train <-  predict(mod,newx=as.matrix(dta$x),type="response")
	    proba.train <-  matrix(c(1-proba.train,proba.train),ncol=2,byrow=FALSE)
	    if (is.null(test)) {proba.test=NULL;predict.test=NULL} 
	    else {
	    proba.test <-  predict(mod,newx=as.matrix(test$x),type="response")
	    proba.test <-  matrix(c(1-proba.test, proba.test),ncol=2,byrow=FALSE)
	    predict.test <-  apply(proba.test,1,which.max)
	    }
	    selected <-  which(abs(coef(mod)[-1])>1e-06) 
	    }
	  else { 
	  	cvmod <-  cv.glmnet(x=as.matrix(dta$x),y=dta$y,family="multinomial",type.measure="class",type.multinomial=type.multinomial,nfolds=n,grouped=FALSE)
	    lambda.min <-  cvmod$lambda.min
	    mod <-  glmnet(x=as.matrix(dta$x),y=dta$y,family="multinomial",lambda=lambda.min,type.multinomial=type.multinomial)
	    proba.train <-  predict(mod,newx=as.matrix(dta$x),type="response")[,,1]
	    if (is.null(test)) {proba.test<- NULL;predict.test<- NULL} 
	    else {
	    proba.test <-  predict(mod,newx=as.matrix(test$x),type="response")[,,1]
	    if (nrow(test$x) == 1) {proba.test <-  matrix(proba.test,nrow=1)}
	    predict.test <-  apply(proba.test,1,which.max)
	    }
	    selected <-  lapply(coef(mod),function(x) which(abs(x[-1])>1e-06 ) ) 
	    } ; 	    return(list(proba.train=proba.train,proba.test=proba.test,predict.test=predict.test,selected=selected,model=mod))
}
	nbclass <- length(unique(faobject$groups))
	
	coefpx <-  LassoML(dta=list(x=faobject$dta,y=faobject$groups),test=test) 
  	proba.test <-  coefpx$proba.test
	
	sB <-  scale(t(faobject $B), center = FALSE, scale = sqrt(faobject$Psi))
    G  <-  solve(diag(ncol(faobject$B)) + tcrossprod(sB))
    sB <-  scale(t(faobject$B), center = FALSE, scale = faobject$Psi)
	
	testzclass <-  lapply(1:nbclass,function(i){tcrossprod(tcrossprod(scale(test$x,center= faobject$meanclass[,i],scale=FALSE),sB),G)})
    testzclass <-  lapply(1:nbclass,function(i){sweep(testzclass[[i]],1,proba.test[,i],"*")})
    testz <-  Reduce('+', testzclass)
    
    fatest <- test$x-tcrossprod(testz, faobject$B) 
	return(list(meanclass = faobject$meanclass, fadta = faobject$fadta, fatest = fatest, 
            Psi = faobject$Psi, B = faobject$B, groups = faobject$groups, dta = NULL))
}
