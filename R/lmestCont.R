lmestCont <- function(responsesFormula = NULL, latentFormula = NULL,
                      data, index, k = 1:4, start = 0,
                      modSel = c("BIC", "AIC"), modBasic = 0,
                      paramLatent = c("multilogit", "difflogit"),
                      weights = NULL, tol = 10^-10,
                      maxit = 5000, out_se = FALSE, output = FALSE,
                      parInit = list(piv = NULL, Pi = NULL,
                                     Mu = NULL, Si = NULL,
                                     Be = NULL, Ga = NULL),
                      fort = TRUE, seed = NULL)
{


  if(inherits(data, "lmestData"))
  {
    data <- data$data
  }else if(!is.data.frame(data))
  {
    data <- as.data.frame(data)
    stop("A data.frame must be provided")
  }

  if(start == 2)
  {
    if(is.null(parInit))
    {
      stop("With start = 2, initial parameters must be provided")

    }else if(!is.null(parInit$Mu))
    {
      k <- dim(parInit$Mu)[2]

    }
  }

  if(!is.null(seed))
  {
    set.seed(seed)
  }

  k <- sort(unique(k))
  nkv <- length(k)

  if(length(index) !=2)
  {
    stop("id and time must be provided")
  }
  id.which <- which(names(data) == index[1])
  tv.which <- which(names(data) == index[2])

  if(length(id.which) == 0)
  {
    stop("the id column does not exist")
  }

  if(length(tv.which) == 0)
  {
    stop("the time column does not exist")
  }

  modSel <- match.arg(modSel, choices = eval(formals(lmestCont)$modSel))
  paramLatent <- match.arg(paramLatent, choices = eval(formals(lmestCont)$paramLatent))

  id <- data[,id.which]
  tv <- data[,tv.which]

if(is.character(id) | is.factor(id))
{
  #warning("id column must be numeric. Coerced in numeric.", call. = FALSE)
  id <- as.numeric(factor(id))
}

if(is.character(tv) | is.factor(tv))
{
  #warning("time column must be numeric. Coerced in numeric.", call. = FALSE)
  tv <- as.numeric(factor(tv))
}

  data.new <- data[,-c(id.which,tv.which), drop = FALSE]
  if(is.null(responsesFormula))
  {
    Y <- data.new
    Xmanifest <- NULL
    Xinitial <- NULL
    Xtrans <- NULL
  }else{
    temp <- getResponses(data = data.new,formula = responsesFormula)
    Y <- temp$Y
    Xmanifest <- temp$X
    Xinitial <- NULL
    Xtrans <- NULL

    Y_names <- colnames(Y)

  }

  if(!is.null(latentFormula))
  {
    temp <- getLatent(data = data.new,latent = latentFormula, responses = responsesFormula)
    Xinitial <- temp$Xinitial
    Xtrans <- temp$Xtrans

    Xinitial_names <- colnames(Xinitial)
    Xtrans_names <- colnames(Xtrans)
  }

  tmp <- long2matrices.internal(Y = Y, id = id, time = tv,
                                Xinitial = Xinitial, Xmanifest = Xmanifest, Xtrans = Xtrans, cont = TRUE)
  model <- tmp$model
  Xinitial <- tmp$Xinitial
  Xtrans <- tmp$Xtrans
  Y <- tmp$Y

  if(is.null(weights))
  {
    freq = tmp$freq
  }else{
    freq = weights
    if(nrow(Y)!=length(weights)) stop("dimensions mismatch between data and weights")
  }

  out = vector("list",nkv)

  dimnames(Y)[[3]] <- Y_names

  if(!is.null(Xinitial))
  {
    dimnames(Xinitial)[[2]] <- Xinitial_names

    if(any(is.na(Xinitial)))
    {
      stop("missing data in the covariates affecting the initial probabilities are not allowed")
    }
  }

  if(!is.null(Xtrans))
  {
    dimnames(Xtrans)[[3]] <- Xtrans_names

    if(any(is.na(Xtrans)))
    {
      stop("missing data in the covariates affecting the transition probabilities are not allowed")
    }
  }


  miss = any(is.na(Y))
  if(miss)
    {
      Yv = cbind(1,Y)
      pYv = prelim.mix(Yv,1)
      thhat = em.mix(prelim.mix(Yv,1))
      rngseed(1)
      Yv = as.matrix(imp.mix(pYv, da.mix(pYv,thhat,steps=100), Yv)[,-1])
      Y = array(Yv,dim(Y))
      cat("Missing data in the dataset. imp.mix function (mix package) used for imputation.\n")
    }



  aicv = rep(NA,nkv)
  bicv = rep(NA,nkv)
  lkv = rep(NA,nkv)


  for(kv in 1:nkv){
    out[[kv]] <- switch(model,
                        "LMbasiccont" = (lmbasic.cont(Y = Y,k = k[kv],start = start,modBasic = modBasic,tol = tol,maxit = maxit,out_se = out_se,
                                                      piv = parInit$piv,Pi = parInit$Pi,Mu = parInit$Mu,Si = parInit$Si, miss = miss, fort = fort)),
                        "LMlatentcont" = (lmcovlatent.cont(Y = Y,X1 = Xinitial,X2 = Xtrans,yv = freq,k = k[kv],
                                                           start = start,tol = tol,maxit = maxit,paramLatent = paramLatent,
                                                           Mu = parInit$Mu, Si = parInit$Si,Be = parInit$Be,Ga = parInit$Ga,output = output, out_se = out_se, miss = miss)))




    aicv[kv] = out[[kv]]$aic
    bicv[kv] = out[[kv]]$bic
    lkv[kv] = out[[kv]]$lk
  }

  if(modSel == "BIC")
  {
    best <- out[[which.min(bicv)]]
  }else if(modSel == "AIC"){
    best <- out[[which.min(aicv)]]
  }

  Bic <- bicv
  Aic <- aicv
  Lk <- lkv

  names(Bic) <- paste("k",k,sep = "=")
  names(Aic) <- paste("k",k,sep = "=")
  names(Lk) <- paste("k",k,sep = "=")

  # if(nkv > 1)
  # {
  #   Bic = bicv
  #   Aic = aicv
  # }
  out <- do.call(c, list(best,
                         list(Lk = lkv, Bic = Bic, Aic = Aic, call = match.call(),data = data)))
  attributes(out)$responsesFormula = responsesFormula
  attributes(out)$latentFormula = latentFormula
  attributes(out)$whichid = id.which
  attributes(out)$whichtv = tv.which
  attributes(out)$id = id
  attributes(out)$time = tv

  class(out) <- class(best)
  return(out)
}
