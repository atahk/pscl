library("devtools")
library("revdepcheck")
library("crancache")
library("parallel")

revdep_reset()
revdep.env <- revdep_env_vars()
revdep.env["R_COMPILE_AND_INSTALL_PACKAGES"] <- "always"
revdep_check(quiet=TRUE, bioc=FALSE,
             timeout = as.difftime(90, units = "mins"),
             num_workers = 4,
             env = revdep.env)
revdep_report()
revdep_report_cran()
