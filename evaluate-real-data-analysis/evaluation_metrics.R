if (outcome == "continuous") {
  evaluation = function(dat){
    model = lm(y ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + age + sex, data = dat)
    dat = dat %>% mutate(residual = model$residuals)
    model = lm(residual ~ PRS, data = dat)
    return(summary(model)$r.squared)
  }
}
if (outcome == "disease") {
  library(RISCA)
  evaluation = function(dat){
    roc_obj = roc.binary(status = "y", variable = "PRS", confounders = ~age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                         data = dat, precision = seq(0.05,0.95, by=0.01))
    return(2*(qnorm(roc_obj$auc)**2))
  }
}