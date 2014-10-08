
FADA = function (faobject, K=10,B=20, nbf.cv = NULL, method = c("glmnet", 
    "sda", "sparseLDA"), sda.method = c("lfdr", "HC"), stop.par = 10, 
    lambda, lambda.var, lambda.freqs, diagonal = FALSE, alpha = 0.1, 
    nfolds = 10) 
{
    nbclass <- length(unique(faobject$groups))
    method <- match.arg(method)
    sda.method <- match.arg(sda.method)
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
       if (nbclass == 2) { proba.train <- matrix(c(1 - proba.train, proba.train), ncol = 2, byrow = FALSE) ; selected <- which(abs(coef(mod)[-1]) > 1e-06)} 
       if (nbclass > 2) {proba.train <- proba.train[, , 1] ; selected <- lapply(coef(mod), function(x) which(abs(x[-1]) > 1e-06))}
        return(list(proba.train = proba.train,model = mod,selected=selected))
    }
    FADA.tmp <- function(faobject, method, sda.method, stop.par, 
        lambda, lambda.var, lambda.freqs, diagonal, alpha) {
        fadta <- faobject$fa.training
        fatest <- faobject$fa.testing
        groups <- faobject$groups
        p <- ncol(faobject$fa.training)
        if (method == "glmnet") {
            out <- LassoML(list(x = fadta, y = groups))
            selected <- out$selected
            proba.test <- predict(out$mod,newx=as.matrix(fatest),type="response")
            if (nbclass == 2) {proba.test <- matrix(c(1 - proba.test, proba.test), 
                ncol = 2, byrow = FALSE)}
            predict.test <- apply(proba.test, 1, which.max)
            out <- out$model
            proba.train <- predict(out,fadta,type="response")
        }
        if (method == "sda") {
            ranking.LDA <- sda::sda.ranking(fadta, groups, diagonal = diagonal, 
                verbose = FALSE)
            if (sda.method == "lfdr") {
                selected <- as.numeric(ranking.LDA[ranking.LDA[, 
                  "lfdr"] < 0.8, "idx"])
            }
            else {
                thr <- which.max(ranking.LDA[1:round(alpha * 
                  p), "HC"])
                selected <- as.numeric(ranking.LDA[1:thr, "idx"])
            }
            out <- sda::sda(fadta[, selected,drop=FALSE], groups, verbose = FALSE)
            pred <- sda::predict.sda(out, fatest[, selected, 
                drop = FALSE], verbose = FALSE)
            proba.test <- pred$posterior
            predict.test <- pred$class
            proba.train <- sda::predict.sda(out,fadta[,selected,drop=FALSE],verbose=FALSE)$posterior
        }
        if (method == "sparseLDA") {
            Xc <- normalize(fadta)
            Xn <- Xc$Xc
            out <- sparseLDA::sda(Xn, factor(groups), stop = -stop.par)
            Xctest <- normalizetest(fatest, Xc)
            Xctest <- matrix(Xctest, nrow = nrow(fatest), byrow = FALSE)
            colnames(Xctest) <- colnames(Xn)
            pred <- sparseLDA::predict.sda(out, Xctest)
            selected <- out$varIndex
            proba.test <- pred$posterior
            predict.test <- pred$class
            proba.train <- sparseLDA::predict.sda(out,Xn)$posterior
        }
        return(list(method = method, selected = selected, proba.train=proba.train,proba.test = proba.test, 
            predict.test = predict.test, mod=out))
    } 
   cv.FADA <- function(train.x, train.y, test.x, test.y,nbf.cv,method = method, sda.method = sda.method, 
                stop.par = stop.par, lambda = lambda, lambda.var = lambda.var, 
                lambda.freqs = lambda.freqs, diagonal = diagonal, 
                alpha = alpha){
   	 fa.train <- decorrelate.train(list(x=train.x,y=train.y),nbf=nbf.cv,verbose=FALSE)
   	 fa.test <- decorrelate.test(fa.train,list(x=test.x))
   	 fada <- FADA.tmp(fa.test, method, sda.method, stop.par, 
        lambda, lambda.var, lambda.freqs, diagonal, alpha)
     return(mean(fada$predict.test != test.y))
   }
    if (! is.null(faobject$fa.testing)) {
        out <- FADA.tmp(faobject, method, sda.method, stop.par, 
            lambda, lambda.var, lambda.freqs, diagonal, alpha)
        proba.train <- out$proba.train
        proba.test <- out$proba.test
        predict.test <- out$predict.test
        selected <- out$selected
        out <- out$mod
        cv.error <- NULL
        cv.error.se <- NULL
    }
    else {
        fadta <- faobject$fa.training
        groups <- faobject$groups
        p <- ncol(fadta)
        if (method == "glmnet") {
            out <- LassoML(list(x = fadta, y = groups))
            selected <- out$selected
           out <- out$model
           proba.train <- predict(out,fadta,type="response")
        }
        if (method == "sda") {
            ranking.LDA <- sda::sda.ranking(fadta, groups, diagonal = diagonal, 
                verbose = FALSE)
            if (sda.method == "lfdr") {selected <- as.numeric(ranking.LDA[ranking.LDA[, "lfdr"] < 0.8, "idx"])}
            if (sda.method == "HC") { thr <- which.max(ranking.LDA[1:round(alpha * p), "HC"]) ;
                selected <- as.numeric(ranking.LDA[1:thr, "idx"]) }
            out <- sda::sda(fadta[, selected,drop=FALSE], groups, verbose = FALSE)
            proba.train <- sda::predict.sda(out,fadta[,selected,drop=FALSE],verbose=FALSE)$posterior
        }
        if (method == "sparseLDA") {
            Xc <- normalize(fadta)
            Xn <- Xc$Xc
            out <- sparseLDA::sda(Xn, factor(groups), stop = -stop.par)
            selected <- out$varIndex
            proba.train <- sparseLDA::predict.sda(out,Xn)$posterior
        }
       cv.out <- crossval(cv.FADA, faobject$data.train,faobject$groups,nbf.cv=nbf.cv,method = method, sda.method = sda.method, 
                stop.par = stop.par, lambda = lambda, lambda.var = lambda.var, 
                lambda.freqs = lambda.freqs, diagonal = diagonal, 
                alpha = alpha,K=K,B=B)
        proba.test <- NULL
        predict.test <- NULL
        cv.error <- cv.out$stat
        cv.error.se <- cv.out$stat.se
    }    
    return(list(method = method, selected = selected, proba.train=proba.train,proba.test = proba.test, 
        predict.test = predict.test, cv.error = cv.error,cv.error.se=cv.error.se, mod=out))
}
