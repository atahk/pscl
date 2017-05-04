odTest <- function(glmobj,
                   alpha=.05,
                   digits=max(3,getOption("digits")-3))
{
  if(class(glmobj)[1]!="negbin")
    stop("this function only works for objects of class negbin\n")

  if(alpha>1 | alpha<0)
    stop("invalid value for alpha\n")
  
  poissonGLM <- glm(formula=eval(glmobj$call$formula),
                    data=eval(glmobj$call$data),
                    family="poisson")
  ## require(stats)
  llhPoisson <- logLik(poissonGLM)
  llhNB <- logLik(glmobj)
  
  d <- 2*(llhNB - llhPoisson)

  ## n.b., distribution of test-statistics is non-standard
  ## see Cameron and Trivedi 1998 p78
  critval <- qchisq(1-(2*alpha), df = 1)
  pval <- pchisq(d, df = 1, lower.tail=FALSE)/2
  
  cat("Likelihood ratio test of H0: Poisson, as restricted NB model:\n")
  cat("n.b., the distribution of the test-statistic under H0 is non-standard\n")
  cat("e.g., see help(odTest) for details/references\n\n")
  
  cat(paste("Critical value of test statistic at the alpha=",
            round(alpha,digits),
            "level:",
            round(critval,digits),
            "\n"))
  cat(paste("Chi-Square Test Statistic = ",
            round(d,digits),
            "p-value =",
            format.pval(pval,digits=digits),
            "\n"))
  
  invisible(NULL)
}
