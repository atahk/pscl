### R code from vignette source 'countreg.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("sandwich")
library("lmtest")
library("MASS")
library("car")
library("pscl")
load("DebTrivedi.rda")

clog <- function(x) log(x + 0.5)
cfac <- function(x, breaks = NULL) {
  if(is.null(breaks)) breaks <- unique(quantile(x, 0:10/10))
  x <- cut(x, breaks, include.lowest = TRUE, right = FALSE)
  levels(x) <- paste(breaks[-length(breaks)], ifelse(diff(breaks) > 1,
    c(paste("-", breaks[-c(1, length(breaks))] - 1, sep = ""), "+"), ""), sep = "")
  return(x)
}

options(prompt = "R> ")
refit_models <- TRUE


###################################################
### code chunk number 2: dt
###################################################
dt <- DebTrivedi[, c(1, 6:8, 13, 15, 18)]


###################################################
### code chunk number 3: dt2
###################################################
dt2 <- DebTrivedi[, -(2:6)]
dt2$region <- relevel(dt2$region, "other")


###################################################
### code chunk number 4: ofp-plot (eval = FALSE)
###################################################
## plot(table(dt$ofp))


###################################################
### code chunk number 5: ofp-plot2
###################################################
plot(table(dt$ofp), xlab = "Number of physician office visits", ylab = "Frequency", axes = FALSE)
axis(2)
axis(1, at = 0:18 * 5, labels = FALSE)
axis(1, at = 0:9 * 10)


###################################################
### code chunk number 6: bad-good
###################################################
par(mfrow = c(1, 2))
plot(ofp ~ numchron, data = dt)
plot(clog(ofp) ~ cfac(numchron), data = dt)


###################################################
### code chunk number 7: bad (eval = FALSE)
###################################################
## plot(ofp ~ numchron, data = dt)


###################################################
### code chunk number 8: clog
###################################################
clog <- function(x) log(x + 0.5)


###################################################
### code chunk number 9: cfac
###################################################
cfac <- function(x, breaks = NULL) {
  if(is.null(breaks)) breaks <- unique(quantile(x, 0:10/10))
  x <- cut(x, breaks, include.lowest = TRUE, right = FALSE)
  levels(x) <- paste(breaks[-length(breaks)], ifelse(diff(breaks) > 1,
    c(paste("-", breaks[-c(1, length(breaks))] - 1, sep = ""), "+"), ""),
    sep = "")
  return(x)
}


###################################################
### code chunk number 10: good (eval = FALSE)
###################################################
## plot(clog(ofp) ~ cfac(numchron), data = dt)


###################################################
### code chunk number 11: ofp2-plot1
###################################################
par(mfrow = c(3, 2))
plot(clog(ofp) ~ health, data = dt, varwidth = TRUE,
  ylab = "Physician office visits (in clogs)", xlab = "Self-perceived health status", main = "health")
plot(clog(ofp) ~ cfac(numchron), data = dt,
  ylab = "Physician office visits (in clogs)", xlab = "Number of chronic conditions", main = "numchron")
plot(clog(ofp) ~ privins, data = dt, varwidth = TRUE,
  ylab = "Physician office visits (in clogs)", xlab = "Covered by private insurance", main = "privins")
plot(clog(ofp) ~ cfac(hosp, c(0:2, 8)), data = dt,
  ylab = "Physician office visits (in clogs)", xlab = "Number of hospital stays", main = "hosp")
plot(clog(ofp) ~ gender, data = dt, varwidth = TRUE,
  ylab = "Physician office visits (in clogs)", xlab = "Gender", main = "gender")
plot(cfac(ofp, c(0:2, 4, 6, 10, 100)) ~ school, data = dt, breaks = 9,
  ylab = "Physician office visits (number of visits)", xlab = "Number of years of education", main = "school")


###################################################
### code chunk number 12: ofp2 (eval = FALSE)
###################################################
## plot(clog(ofp) ~ health, data = dt, varwidth = TRUE)
## plot(clog(ofp) ~ cfac(numchron), data = dt)
## plot(clog(ofp) ~ privins, data = dt, varwidth = TRUE)
## plot(clog(ofp) ~ cfac(hosp, c(0:2, 8)), data = dt)
## plot(clog(ofp) ~ gender, data = dt, varwidth = TRUE)
## plot(cfac(ofp, c(0:2, 4, 6, 10, 100)) ~ school, data = dt, breaks = 9)


###################################################
### code chunk number 13: models
###################################################
if(refit_models & file.exists("countreg-models.rda")) file.remove("countreg-models.rda")
if(file.exists("countreg-models.rda")) {
  load("countreg-models.rda")
} else {
  fm_pois   <-      glm(ofp ~ ., data = dt, family = poisson)
  fm_qpois  <-      glm(ofp ~ ., data = dt, family = quasipoisson)
  fm_nbin   <-   glm.nb(ofp ~ ., data = dt)
  fm_zinb0  <- zeroinfl(ofp ~ ., data = dt, dist = "negbin")
  fm_zinb   <- zeroinfl(ofp ~ . | hosp + numchron + privins + school + gender, data = dt, dist = "negbin")
  fm_hurdle0<-   hurdle(ofp ~ ., data = dt, dist = "negbin")
  fm_hurdle <-   hurdle(ofp ~ . | hosp + numchron + privins + school + gender, data = dt, dist = "negbin")
  fm_hurdle2<-   hurdle(ofp ~ ., data = dt2, dist = "negbin")
  if(!refit_models) save(fm_pois, fm_qpois, fm_nbin, fm_zinb0, fm_zinb, fm_hurdle0, fm_hurdle, fm_hurdle2, file = "countreg-models.rda")
}


###################################################
### code chunk number 14: poisson (eval = FALSE)
###################################################
## fm_pois <- glm(ofp ~ ., data = dt, family = poisson)


###################################################
### code chunk number 15: summary-poisson
###################################################
summary(fm_pois)


###################################################
### code chunk number 16: coeftest-poisson
###################################################
coeftest(fm_pois, vcov = sandwich)


###################################################
### code chunk number 17: quasipoisson (eval = FALSE)
###################################################
## fm_qpois <- glm(ofp ~ ., data = dt, family = quasipoisson)


###################################################
### code chunk number 18: summary-quasipoisson (eval = FALSE)
###################################################
## summary(fm_qpois)


###################################################
### code chunk number 19: nbin (eval = FALSE)
###################################################
## fm_nbin <- glm.nb(ofp ~ ., data = dt)
## summary(fm_nbin)


###################################################
### code chunk number 20: hurdle0 (eval = FALSE)
###################################################
## fm_hurdle0 <- hurdle(ofp ~ ., data = dt, dist = "negbin")


###################################################
### code chunk number 21: summary-hurdle0
###################################################
summary(fm_hurdle0)


###################################################
### code chunk number 22: hurdle (eval = FALSE)
###################################################
## fm_hurdle <- hurdle(ofp ~ . | hosp + numchron + privins + school + gender,
##   data = dt, dist = "negbin")


###################################################
### code chunk number 23: waldtest-hurdle
###################################################
waldtest(fm_hurdle0, fm_hurdle)


###################################################
### code chunk number 24: lrtest-hurdle (eval = FALSE)
###################################################
## lrtest(fm_hurdle0, fm_hurdle)


###################################################
### code chunk number 25: summary-table
###################################################
fm <- list("ML-Pois" = fm_pois, "Adj-Pois" = fm_pois, "Quasi-Pois" = fm_qpois, "NB" = fm_nbin,
  "Hurdle-NB" = fm_hurdle, "ZINB" = fm_zinb)
fm_summary <- matrix(character(6 * 33), ncol = 6)
colnames(fm_summary) <- names(fm)
rownames(fm_summary) <- c(as.vector(rbind(names(coef(fm_hurdle, model = "count")), "")),
  as.vector(rbind(names(coef(fm_hurdle, model = "zero")), "")),
  "no.\\ parameters", "$\\log L$", "AIC", "BIC", "$\\sum_i \\hat f_i(0)$")
rownames(fm_summary)[1:28] <- ifelse(rownames(fm_summary)[1:28] == "", "", 
  paste("\\code{", rownames(fm_summary)[1:28], "}", sep = ""))
fm_summary[1:8 * 2 - 1,] <- sapply(fm,
  function(x) paste("$", format(round(coef(x)[1:8], digits = 3)), "$\\phantom{)}", sep = ""))
fm_summary[1:8 * 2,] <- sapply(
  c(list("ML-Pois" = vcov(fm_pois), "Adj-Pois" = sandwich(fm_pois)),
  lapply(fm[-(1:2)], function(x) vcov(x))), function(x)
  paste("(", format(round(sqrt(diag(x))[1:8], digits = 3)), ")", sep = ""))
fm_summary[1:6 * 2 + 15,] <- cbind(NA, NA, NA, NA, sapply(fm[5:6], function(x)
  paste("$", format(round(coef(x, model = "zero"), digits = 3)), "$\\phantom{)}", sep = "")))
fm_summary[1:6 * 2 + 16,] <- cbind(NA, NA, NA, NA, sapply(fm[5:6],
  function(x) paste("(", format(round(sqrt(diag(vcov(x)))[-(1:8)], digits = 3)),  ")", sep = "")))
fm_summary[29,] <- sapply(fm, function(x) attr(logLik(x), "df"))
fm_summary[30,] <- paste("$", format(sapply(fm, function(x) round(logLik(x), digits = 1))), "$", sep = "")
fm_summary[31,] <- format(round(sapply(fm, AIC), digits = 1))
fm_summary[32,] <- format(round(sapply(fm, AIC, k = log(nrow(dt))), digits = 1))
fm_summary[33,] <- round(c("ML-Pois" = sum(dpois(0, fitted(fm_pois))),
  "Adj-Pois" = NA,
  "Quasi-Pois" = NA,
  "NB" = sum(dnbinom(0, mu = fitted(fm_nbin), size = fm_nbin$theta)),
  "NB-Hurdle" = sum(predict(fm_hurdle, type = "prob")[,1]),
  "ZINB" = sum(predict(fm_zinb, type = "prob")[,1])))
fm_summary[30:33,2:3] <- NA
fm_summary[is.na(fm_summary)] <- " "
fm_summary <- paste(apply(cbind(rownames(fm_summary), fm_summary), 1, paste, collapse = " & "), "\\\\")
fm_summary[c(16, 28, 33)] <- paste(fm_summary[c(16, 28, 33)], "\\hline")
writeLines(fm_summary)


###################################################
### code chunk number 26: zinb0 (eval = FALSE)
###################################################
## fm_zinb0 <- zeroinfl(ofp ~ ., data = dt, dist = "negbin")


###################################################
### code chunk number 27: zinb (eval = FALSE)
###################################################
## fm_zinb <- zeroinfl(ofp ~ . | hosp + numchron + privins + school + gender,
##   data = dt, dist = "negbin")


###################################################
### code chunk number 28: waldtest-zinb
###################################################
waldtest(fm_zinb0, fm_zinb)


###################################################
### code chunk number 29: summary-zinb (eval = FALSE)
###################################################
## summary(fm_zinb)


###################################################
### code chunk number 30: coef-count
###################################################
fm <- list("ML-Pois" = fm_pois, "Quasi-Pois" = fm_qpois, "NB" = fm_nbin,
  "Hurdle-NB" = fm_hurdle, "ZINB" = fm_zinb)
sapply(fm, function(x) coef(x)[1:8])


###################################################
### code chunk number 31: se-count
###################################################
cbind("ML-Pois" = sqrt(diag(vcov(fm_pois))),
  "Adj-Pois" = sqrt(diag(sandwich(fm_pois))),
  sapply(fm[-1], function(x) sqrt(diag(vcov(x)))[1:8]))


###################################################
### code chunk number 32: logLik
###################################################
rbind(logLik = sapply(fm, function(x) round(logLik(x), digits = 0)),
  Df = sapply(fm, function(x) attr(logLik(x), "df")))


###################################################
### code chunk number 33: zero-counts
###################################################
round(c("Obs" = sum(dt$ofp < 1),
  "ML-Pois" = sum(dpois(0, fitted(fm_pois))),
  "NB" = sum(dnbinom(0, mu = fitted(fm_nbin), size = fm_nbin$theta)),
  "NB-Hurdle" = sum(predict(fm_hurdle, type = "prob")[,1]),
  "ZINB" = sum(predict(fm_zinb, type = "prob")[,1])))


###################################################
### code chunk number 34: coef-zero
###################################################
t(sapply(fm[4:5], function(x) round(x$coefficients$zero, digits = 3)))


###################################################
### code chunk number 35: dt2a (eval = FALSE)
###################################################
## dt2 <- DebTrivedi[, -(2:6)]
## dt2$region <- relevel(dt2$region, "other")


###################################################
### code chunk number 36: hurdle2 (eval = FALSE)
###################################################
## fm_hurdle2 <- hurdle(ofp ~ ., data = dt2, dist = "negbin")


###################################################
### code chunk number 37: hurdle2-summary
###################################################
cfz <- coef(fm_hurdle2, model = "zero")
cfc <- coef(fm_hurdle2, model = "count")
se <- sqrt(diag(sandwich(fm_hurdle2)))
round(cbind(zero = cfz, zero_t = cfz/se[-seq(along = cfc)], 
  count = cfc, count_t = cfc/se[seq(along = cfc)]),
  digits = 3)[c(3, 2, 4, 5, 7, 6, 8, 9:17, 1),]
logLik(fm_hurdle2)
1/fm_hurdle2$theta


