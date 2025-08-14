######################################################################
############## feature selection
library(randomForestSRC)
library(CoxBoost)
library(MASS)
# Lasso --------------------------------------------------------------
Lasso.fs <- function(train) {
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      nfolds = 10,
                      family = "cox", # cox
                      grouped = TRUE, 
                      nlambda=100,
                      alpha = 1, 
                      # type.measure = "mse"
  )
  coef.min = coef(cv.fit, s = "lambda.min") 
  active.min = which(as.numeric(coef.min)!=0)
  rid <- colnames(x2)[active.min]
  return(rid)
}
# RSF --------------------------------------------------------------
RSF.fs <- function(train) {
  mod.RSF <- rfsrc(Surv(time, status) ~ .,
                   data = train,
                   ntree = 1000,
                   splitrule = "logrank",
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = 123456
  )
  vi <- data.frame(imp=vimp.rfsrc(mod.RSF)$importance)
  vi$imp <- (vi$imp-min(vi$imp))/(max(vi$imp)-min(vi$imp))
  vi$ID <- rownames(vi)
  rid <- rownames(vi)[vi$imp>0.4] ## 0.4 can be adjust
  return(rid)
}
# CoxBoost --------------------------------------------------------------
CoxBoost.fs <- function(train) {
  time <- train$time
  status <- train$status
  x <- train %>%
    dplyr::select(-c(1:2)) %>%
    as.matrix()
  # determine penalty parameter
  set.seed(123456)
  optim.CoxBoost <- optimCoxBoostPenalty(
    time = time, status = status, x = x,
    trace = T, start.penalty = 500
  )
  # Fit with obtained penalty parameter and optimal number of boosting
  # steps obtained by cross-validation
  mod.CoxBoost <- CoxBoost(
    time = time, status = status, x = x,
    stepno = optim.CoxBoost$cv.res$optimal.step,
    penalty = optim.CoxBoost$penalty
  )
  rid <- names(coef(mod.CoxBoost)[which(coef(mod.CoxBoost)!=0)])
  return(rid)
}
# stepwiseCox.both --------------------------------------------------------------
stepwiseCox.both.fs <- function(train) {
  fit0 <- coxph(Surv(time, status) ~ ., data = train)
  fit_both <- MASS::stepAIC(fit0, direction = "both")
  rid <- names(coef(fit_both))
  return(rid)
}
