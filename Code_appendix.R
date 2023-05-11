################################################################################
### Code Appendix - RCT Generalization: PARADIGM-HF
##
## Created on: 2023-05-11
## Author: HoJin Shin
################################################################################



### Prepare environment
################################################################################

### Load libraries
library(readr) # To load .csv data.
library(dplyr)
library(rms)
library(glmnet)

### Load data
trial_data <- read_csv("analytic_file_trial.csv", show_col_types = F)
target_data <- read_csv("analytic_file_target.csv", show_col_types = F)

### Calibration function utilizing rms package's calibrate()
calibration.func <- function(fitted.model, survtime, m, b, lower.lim.x,  upper.lim.x, lower.lim.y, upper.lim.y){
  
  model_info <- deparse(substitute(fitted.model))
  
  survival <- survest(fitted.model, times = survtime, conf.int=FALSE)$surv # Predicted survival probabilities
  output1 <- summary(survival)
  output2 <- unique(survival)
  
  ny <- dim(fitted.model$y)
  m <- m
  g <- max(1, floor(ny[1]/m))
  cuts <- unique(quantile(c(0, 1, survival), seq(0, 1, length=g+1), na.rm=TRUE))
  
  set.seed(1119)
  cal.km <- calibrate(fitted.model, cmethod = 'KM', u=survtime, cuts=cuts, B=b)
  output3 <- summary(cal.km)
  cal.hare <- calibrate(fitted.model, cmethod = 'hare', u=survtime, B=b)
  output4 <- summary(cal.hare)
  
  plot(cal.km, xlim = c(lower.lim.x, upper.lim.x), ylim = c(lower.lim.y, upper.lim.y))
  plot(cal.hare, add=T)
  
  return(assign(paste("calibration", model_info, survtime, sep="_"), list(output1, output2, output3, output4), envir = .GlobalEnv))

}


### Set common values
## Parameters for validation and calibration
survTime <- 730 # Change the time based on the research question.

method <- "boot"
bootItr <- 300
cmethod <- "KM"
m <- 1000
llx <- 0.8
ulx <- 1
lly <- 0.8
uly <- 1

## Model formula with all candidate predictors (rms package format)
fmla.full <- Surv(outcome_te, outcome_i) ~ catg(race) + poly(age_sc, 3, raw = T) + sex +  
  BHTN + BDM + BAF + BHHF + BMI + BSTROKE + BCHD + BANG + BANA + BCABG + BPCI + 
  ACE + ARB + DIUR + DIGOX + BETAB + MRA + ANTIPLAT + ANTICOAG + STATIN + LIPID + ICD + CRT + GLPDPP4



###
### Model derivation and validation (primary outcome as an example)
################################################################################

### Sacubitril/valsartan
# Use data from those who received sacubitril/valsartan.
# There should be no missing values for the selected predictors.
ex <- trial_data %>% filter(exposure == "sacubitril/valsartan") 
dd.ex <- datadist(ex)
options(datadist = 'dd.ex')

# Scale continuous variables
ex$age_sc <- scale(ex$age, center = T, scale = T) # scale is generic function whose default method centers and/or scales the columns of a numeric matrix.

# Create predictor and outcome vectors
X <- model.matrix(fmla.full, ex)[,-1]
y <- cbind(time = ex$outcome_te, status = ex$outcome_i)

penalty.factor <- c(rep(0, 5), rep(1, ncol(X)-5)) ## penalty of 0 ==> always in the model: this command makes race and linear 'AGE' stay in the model.

# Select predictors, using the relaxed LASSO
set.seed(1119)
fit.coxph.cv <- cv.glmnet(X, y, family = "cox", penalty.factor = penalty.factor, alpha = 1, nfolds = 10, relax = T)
fit.coxph.cv
plot(fit.coxph.cv)

# Get the nonzero coefficients
selectedVar <- rownames(coef(fit.coxph.cv))[as.matrix(coef(fit.coxph.cv, fit.coxph.cv$relaxed$lambda.min)) != 0]
selectedVar

# Write a formula with the selected variables
formula.ex <- Surv(outcome_te, outcome_i) ~ catg(race) + pol(age, 3) + sex + 
  BHTN + BDM + BHHF + BMI + BSTROKE + BCHD + BANG + BPCI + 
  ACE + DIUR + DIGOX + BETAB + MRA + ANTIPLAT + ANTICOAG + STATIN + ICD + CRT + GLPDPP4

# Validate the final model
model.fit.ex <- cph(formula.ex, data = ex, x = T, y = T, surv = T, time.inc = survTime)
set.seed(1119) 
validate(model.fit.ex, method = method, B = bootItr, bw = FALSE, dxy = T, type = "residual")
cal.km <- calibrate(model.fit.ex, u = survTime, cmethod = cmethod, m = m, B = bootItr)
plot(cal.km)
calibration.func(model.fit.ex, survtime = survTime, m = m, b = bootItr, lower.lim.x = llx, upper.lim.x = ulx, lower.lim.y = lly, upper.lim.y = uly)



### Enalapril
# Use data from those who received sacubitril/valsartan.
# There should be no missing values for the selected predictors.
unex <- trial_data %>% filter(exposure == "sacubitril/valsartan") 
dd.unex <- datadist(unex)
options(datadist = 'dd.unex')

# Scale continuous variables
unex$age_sc <- scale(unex$age, center = T, scale = T) # scale is generic function whose default method centers and/or scales the columns of a numeric matrix.

# Create predictor and outcome vectors
X <- model.matrix(fmla.full, unex)[,-1]
y <- cbind(time = unex$outcome_te, status = unex$outcome_i)

penalty.factor <- c(rep(0, 5), rep(1, ncol(X)-5)) ## penalty of 0 ==> always in the model: this command makes race and linear 'AGE' stay in the model.

# Select predictors, using the relaxed LASSO
set.seed(1119)
fit.coxph.cv <- cv.glmnet(X, y, family = "cox", penalty.factor = penalty.factor, alpha = 1, nfolds = 10, relax = T)
fit.coxph.cv
plot(fit.coxph.cv)

# Get the nonzero coefficients
selectedVar <- rownames(coef(fit.coxph.cv))[as.matrix(coef(fit.coxph.cv, fit.coxph.cv$relaxed$lambda.min)) != 0]
selectedVar

# Write a formula with the selected variables
formula.unex <- Surv(outcome_te, outcome_i) ~ catg(race) + pol(age, 2) + sex + 
  BDM + BAF + BHHF + BMI + BANG + BCABG + DIUR + DIGOX + BETAB + LIPID

# Validate the final model
model.fit.unex <- cph(formula.unex, data = unex, x = T, y = T, surv = T, time.inc = survTime)
set.seed(1119) 
validate(model.fit.unex, method = method, B = bootItr, bw = FALSE, dxy = T, type = "residual")
cal.km <- calibrate(model.fit.unex, u = survTime, cmethod = cmethod, m = m, B = bootItr)
plot(cal.km)
calibration.func(model.fit.unex, survtime = survTime, m = m, b = bootItr, lower.lim.x = llx, upper.lim.x = ulx, lower.lim.y = lly, upper.lim.y = uly)



###
### Generating 10,000 bootstrap model fits in the trial data
################################################################################

## Number of Bootstrap resampling
B <- 10000

set.seed(1119)
results = NULL
fit.list <- list()
baseRisk.list <- list()
for (i in 1:B){
  
  ### Sacubitril/valsartan
  
  # Draw a bootstrap sample from the original trial data
  ex.boot <- ex[sample(nrow(ex), nrow(ex), TRUE), ]
  fit.ex.boot <- cph( formula.ex, data = ex.boot, x = T, y = T, surv = T, time.inc = survTime)
  
  basehaz.ex <- basehaz(coxph(formula.ex, data = df.ex.boot), centered = FALSE)
  baseRisk.ex <- basehaz.ex[basehaz.ex$time == 730, ]$hazard
  
  # Survival probabilities
  modmat.ex <- model.matrix(formula.ex, data = df.ex.boot)[, -1]
  eta.ex <- modmat.ex %*% fit.ex.boot$coefficients # Linear predictor = b1*x1 + b2*x2 + ...
  pred.ex <- exp(-baseRisk.ex * exp(eta.ex))
  
  ### Enalapril
  
  # Draw a bootstrap sample from the original trial data
  unex.boot <- unex[sample(nrow(unex), nrow(unex), TRUE), ]
  fit.unex.boot <- cph( formula.unex, data = unex.boot, x = T, y = T, surv = T, time.inc = survTime)
  
  basehaz.unex <- basehaz(coxph(formula.unex, data = df.unex.boot), centered = FALSE)
  baseRisk.unex <- basehaz.unex[basehaz.unex$time == 730, ]$hazard
  
  # Survival probabilities
  modmat.unex <- model.matrix( formula.unex[-2], unex.boot)[,-1]
  eta.unex <- modmat.unex %*% fit.unex.boot$coefficients # linear predictor
  pred.unex <- exp(-baseRisk.unex * exp(eta.unex))
  
  # Sacubitril/valsartan vs. enalapril
  res <- cbind(
    
    risk.ex = mean(1-pred.ex), # Mean predicted probabilities of the outcome for participants treated with sacubitril/valsartan
    risk.unex = mean(1-pred.unex), # Mean predicted probabilities of the outcome for participants treated with enalapril
    rr = mean(1-pred.ex)/mean(1-pred.unex),
    rd = mean(1-pred.ex)-mean(1-pred.unex)
    
  )
  
  results = rbind(results, res)
  fit.list[[i]] <- list(fit.ex.boot$coefficients, fit.unex.boot$coefficients)
  baseRisk.list[[i]] <- list(baseRisk.ex, baseRisk.unex)
  
}

### Obtain bootstrap mean and confidence limits for the trial data
results_prim_2yr <- as.data.frame(results) %>% mutate(
  
  RISK.EX = paste0(format(round(mean(results$risk.ex)*100, 2), nsmall = 2), " (", format(round(quantile(results$risk.ex, c(.025, .975))[1]*100, 2), nsmall = 2), ", ", format(round(quantile(results$risk.ex, c(.025, .975))[2]*100, 2), nsmall = 2), ")"),
  RISK.UNEX = paste0(format(round(mean(results$risk.unex)*100, 2), nsmall = 2), " (", format(round(quantile(results$risk.unex, c(.025, .975))[1]*100, 2), nsmall = 2), ", ", format(round(quantile(results$risk.unex, c(.025, .975))[2]*100, 2), nsmall = 2), ")"),
  RR = paste0(format(round(mean(results$rr), 2), nsmall = 2), " (", format(round(quantile(results$rr, c(.025, .975))[1], 2), nsmall = 2), ", ", format(round(quantile(results$rr, c(.025, .975))[2], 2), nsmall = 2), ")"),
  RD = paste0(format(round(mean(results$rd)*100, 2), nsmall = 2), " (", format(round(quantile(results$rd, c(.025, .975))[1]*100, 2), nsmall = 2), ", ", format(round(quantile(results$rd, c(.025, .975))[2]*100, 2), nsmall = 2), ")")
  
) %>% dplyr::select(RISK.EX:RD) %>% distinct_all


###
### Prediction in the target population
################################################################################

### Estimation function
estimationFunction <- function(outcome, df, B, fit, formula.ex, formula.unex, baseRisk, survTime){
  
  # outcome: name of the outcome
  # df: data
  # B: Number of iteration
  # fit: list of bootstrap models for sacubitril/valsartan and enalapril
  # formula.ex: model formula for sacubitril/valsartan
  # formula.unex: model formula for enalapril
  # baseRisk: list of baseline cumulative hazards at survTime for sacubitril/valsartan and enalapril
  
  outcome_name <- deparse(substitute(outcome))
  
  set.seed(1119)
  results = NULL
  for (i in 1:B){
    
    # Survival probabilities for sacubitril/valsartan
    modmat.ex <- model.matrix( formula.ex[-2], df)[,-1]
    eta.ex <- modmat.ex %*% fit[[i]][[1]] # linear predictor
    pred.ex <- exp(-baseRisk[[i]][[1]] * exp(eta.ex))
    
    # Survival probabilities for enalapril
    modmat.unex <- model.matrix( formula.unex[-2], df)[,-1]
    eta.unex <- modmat.unex %*% fit[[i]][[2]] # linear predictor
    pred.unex <- exp(-baseRisk[[i]][[2]] * exp(eta.unex))
    
    res <- df %>% mutate(
      
      risk.ex = mean(1-pred.ex),
      risk.unex = mean(1-pred.unex),
      rr = risk.ex/risk.unex,
      rd = risk.ex-risk.unex
      
    ) %>% dplyr::select(risk.ex, risk.unex, rr, rd) %>% distinct_all
    
    results = rbind(results, res)
    
  }
  
  res <- results %>% mutate(
    
    RISK.EX = paste0(format(round(mean(results$risk.ex)*100, 2), nsmall = 2), " (", format(round(quantile(results$risk.ex, c(.025, .975))[1]*100, 2), nsmall = 2), ", ", format(round(quantile(results$risk.ex, c(.025, .975))[2]*100, 2), nsmall = 2), ")"),
    RISK.UNEX = paste0(format(round(mean(results$risk.unex)*100, 2), nsmall = 2), " (", format(round(quantile(results$risk.unex, c(.025, .975))[1]*100, 2), nsmall = 2), ", ", format(round(quantile(results$risk.unex, c(.025, .975))[2]*100, 2), nsmall = 2), ")"),
    RR = paste0(format(round(mean(results$rr), 2), nsmall = 2), " (", format(round(quantile(results$rr, c(.025, .975))[1], 2), nsmall = 2), ", ", format(round(quantile(results$rr, c(.025, .975))[2], 2), nsmall = 2), ")"),
    RD = paste0(format(round(mean(results$rd)*100, 2), nsmall = 2), " (", format(round(quantile(results$rd, c(.025, .975))[1]*100, 2), nsmall = 2), ", ", format(round(quantile(results$rd, c(.025, .975))[2]*100, 2), nsmall = 2), ")")
    
  ) %>% dplyr::select(RISK.EX:RD) %>% distinct_all
  
  return(assign(paste(outcome_name, survTime, sep = "."), res, envir = .GlobalEnv))
  
}

estimationFunction(primary, target_data, B, fit.list, formula.ex, formula.unex, baseRisk.list, survTime)





