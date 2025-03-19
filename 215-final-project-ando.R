library(tidyverse)
library(TH.data)
library(survival)
library(survminer)
library(dbplyr)
library(zoo)
library(tidyverse)
library(gtsummary)
library(psych)
library(kableExtra)
library(finalfit)
library(PAmeasures)

data("GBSG2")

data_gbsg2 <- as_tibble(GBSG2)

fit <- survfit(Surv(time, cens) ~ horTh, data = data_gbsg2)
ggsurvplot(fit, data = data_gbsg2, pval = TRUE)

data_gbsg2$tgrade <- factor(data_gbsg2$tgrade, ordered = FALSE)

full.cox.model <- coxph(Surv(time, cens) ~ horTh + tgrade + horTh:tgrade +
                          menostat + age + tsize + pnodes + progrec + estrec,
                        data = data_gbsg2)

# Full model
summary(full.cox.model) 

best.model.aic <- stepAIC(full.cox.model, direction="both", k = 2, trace = 1)

# AIC model
summary(best.model.aic) 

aic.cox.model <- coxph(Surv(time, cens) ~ horTh + tgrade + tsize + pnodes + progrec,
                       data = data_gbsg2, x = TRUE, y = TRUE)

# R^2
pam.coxph(aic.cox.model)

# Likelyhood ratio test
anova(best.model.aic, full.cox.model)

# Concordance
concordance(aic.cox.model)

# Concordance with time weight
concordance(aic.cox.model, timewt="n/G2")

# Weibul AFT
weibreg <- survreg(Surv(time, cens) ~ horTh + tgrade + tsize + pnodes + progrec,
                   data = data_gbsg2, dist = "weibull", x=TRUE, y=TRUE)

# Weibul AFT summary
summary(weibreg)

# Weibul AFT R^2
pam.coxph(weibreg)

# Weibul AFT Concordance
concordance(weibreg)

# Weibul AFT Concordance with time weight
concordance(weibreg, timewt="n/G2")

# Log-normal AFT
logreg <- survreg(Surv(time, cens) ~ horTh + tgrade + tsize + pnodes + progrec,
                  data = data_gbsg2, dist = "lognormal", x=TRUE, y=TRUE)

# Log-normal AFT summary
summary(logreg)

# Log-normal AFT R^2
pam.coxph(logreg)

# Log-normal AFT Concordance
concordance(logreg)

# Log-normal AFT Concordance with time weight
concordance(logreg, timewt="n/G2")

# Log-logistic AFT
loglogreg <- survreg(Surv(time, cens) ~ horTh + tgrade + tsize + pnodes + progrec,
                     data = data_gbsg2, dist = "loglogistic", x=TRUE, y=TRUE)

# Log-logistic AFT summary
summary(loglogreg)

# Log-logistic AFT R^2
pam.coxph(loglogreg)

# Log-logistic AFT Concordance
concordance(loglogreg)

# Log-logistic AFT Concordance with time weight
concordance(loglogreg, timewt="n/G2")

# Exponential AFT
expreg <- survreg(Surv(time, cens) ~ horTh + tgrade + tsize + pnodes + progrec,
                  data = data_gbsg2, dist = "exponential", x=TRUE, y=TRUE)

# Exponential AFT summary
summary(expreg)

# Exponential AFT R^2
pam.coxph(expreg)

# Exponential AFT Concordance
concordance(expreg)

# Exponential AFT Concordance with time weight
concordance(expreg, timewt="n/G2")

# Cox-Snell residuals
mg.residual <- resid(best.model.aic, type ="martingale")

cs.residual <- data_gbsg2$cens - mg.residual

fit.cs <-survfit(Surv(cs.residual, data_gbsg2$cens) ~ 1)

H.cs <- cumsum(fit.cs$n.event / fit.cs$n.risk)

plot(fit.cs$time, H.cs,type = 's',col = 'blue',
     xlab = 'Residual',
     ylab = 'Nelson-Aalen Cum. Hazard',
     xlim = c(0, 2))

abline(0, 1, col = 'red', lty =2)

# Anderson residuals
fit <- basehaz(coxph(Surv(time, cens) ~ horTh + strata(factor(tgrade)) +
                       tsize + pnodes + progrec, data = data_gbsg2,
                     ties = 'breslow'), centered = TRUE)
stage0 <- data.frame(
  "H1" = fit$hazard[fit$strata == "I"],
  "time" = fit$time[fit$strata == "I"])
stage1 <- data.frame(
  "H2" = fit$hazard[fit$strata == "II"],
  "time" = fit$time[fit$strata == "II"])
stage2 <- data.frame(
  "H3" = fit$hazard[fit$strata == "III"],
  "time" = fit$time[fit$strata == "III"])

join1 <- full_join(stage0, stage1, by = "time") %>%
  arrange(time) %>% do(na.locf(.))

join2 <- full_join(join1, stage2, by = "time") %>%
  arrange(time) %>% do(na.locf(.))


plot(join2$H1 ~ join2$H2, 
     type = 's', xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'I', ylab = 'II')

plot(join2$H2 ~ join2$H3, 
     type = 's', xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'II', ylab = 'III')

plot(join2$H1 ~ join2$H3,
     type = 's', xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'I', ylab = 'III')
