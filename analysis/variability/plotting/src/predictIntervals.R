## computes confidence intervals for mixed models
## taktes model, prediction data, name of the response variable, and
## the family of the model so the inverse link function can be applied

predict.int <- function(mod,
                        dd,
                        y,
                        family="gaussian"){
  ## matrix-multiply X by the parameter vector β to get the predictions
  ## (or linear predictor in the case of GLM(M)s);
  mm <- model.matrix(terms(mod), dd)
  dd[,y] <- mm %*% fixef(mod)
  ## extract the variance-covariance matrix of the parameters V
  ## compute XVX′^-1 to get the variance-covariance matrix of the
  ## predictions and extract the diagonal of this matrix to get
  ## variances of predictions
  pvar1 <- diag(mm %*% vcov(mod) %*% t(mm))
  ## tvar1 <- pvar1 + VarCorr(mod)$Plot[1] + VarCorr(mod)$Transect[1] +
  ##   VarCorr(mod)$SurveyYear[1]
  ## take the square-root of the variances to get the standard
  ## deviations (errors) of the predictions;
  ## compute confidence intervals based on a Normal approximation
  new.dd <- data.frame(dd,
                       plo=dd[,y] - 2*sqrt(pvar1),
                       phi=dd[,y] + 2*sqrt(pvar1))
  ## run the confidence interval boundaries through the inverse-link
  ## function.
  if(family=="binomial"){
    new.dd[,y] <- inv.logit(new.dd[,y])
    new.dd$plo <- inv.logit(new.dd$plo)
    new.dd$phi <- inv.logit(new.dd$phi)
  } else if(family=="poisson"){
    new.dd[,y] <- exp(new.dd[,y])
    new.dd$plo <- exp(new.dd$plo)
    new.dd$phi <- exp(new.dd$phi)
   }
  return(new.dd)
}
