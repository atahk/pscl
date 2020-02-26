# AICcmodavg

<details>

* Version: 2.2-2
* Source code: https://github.com/cran/AICcmodavg
* Date/Publication: 2019-05-29 21:20:03 UTC
* Number of recursive dependencies: 80

Run `revdep_details(,"AICcmodavg")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Packages unavailable to check Rd xrefs: ‘multcomp’, ‘gmodels’
    ```

# BSGS

<details>

* Version: 2.0
* Source code: https://github.com/cran/BSGS
* Date/Publication: 2015-06-24 01:07:18
* Number of recursive dependencies: 5

Run `revdep_details(,"BSGS")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    BSGS.Sample: no visible global function definition for ‘rgamma’
    BSGS.Simple: no visible global function definition for ‘rgamma’
    CompWiseGibbs: no visible global function definition for ‘runif’
    CompWiseGibbs: no visible global function definition for ‘rnorm’
    CompWiseGibbsSMP: no visible global function definition for ‘rbinom’
    CompWiseGibbsSMP: no visible global function definition for ‘runif’
    CompWiseGibbsSMP: no visible global function definition for ‘rnorm’
    CompWiseGibbsSimple: no visible global function definition for ‘runif’
    CompWiseGibbsSimple: no visible global function definition for ‘rnorm’
    CompWiseGibbsSimple: no visible global function definition for ‘rgamma’
    GroupMH: no visible global function definition for ‘runif’
    GroupZ.MH: no visible global function definition for ‘runif’
    Undefined global functions or variables:
      rbinom rgamma rnorm runif
    Consider adding
      importFrom("stats", "rbinom", "rgamma", "rnorm", "runif")
    to your NAMESPACE file.
    ```

# congressbr

<details>

* Version: 0.2.2
* Source code: https://github.com/cran/congressbr
* Date/Publication: 2019-12-12 11:20:02 UTC
* Number of recursive dependencies: 110

Run `revdep_details(,"congressbr")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘devtools’
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 1 marked UTF-8 string
    ```

# Countr

<details>

* Version: 3.5.4
* Source code: https://github.com/cran/Countr
* URL: https://github.com/GeoBosh/Countr https://geobosh.github.io/Countr/
* BugReports: https://github.com/GeoBosh/Countr/issues
* Date/Publication: 2019-08-21 11:30:06 UTC
* Number of recursive dependencies: 101

Run `revdep_details(,"Countr")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 14.5Mb
      sub-directories of 1Mb or more:
        doc    2.0Mb
        libs  11.6Mb
    ```

# DClusterm

<details>

* Version: 1.0-1
* Source code: https://github.com/cran/DClusterm
* Date/Publication: 2020-02-25 13:10:06 UTC
* Number of recursive dependencies: 57

Run `revdep_details(,"DClusterm")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘INLA’
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 2 marked Latin-1 strings
    ```

# emIRT

<details>

* Version: 0.0.11
* Source code: https://github.com/cran/emIRT
* Date/Publication: 2020-02-04 06:20:02 UTC
* Number of recursive dependencies: 4

Run `revdep_details(,"emIRT")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 42.1Mb
      sub-directories of 1Mb or more:
        data   2.0Mb
        libs  39.9Mb
    ```

# glmmTMB

<details>

* Version: 1.0.0
* Source code: https://github.com/cran/glmmTMB
* URL: https://github.com/glmmTMB
* BugReports: https://github.com/glmmTMB/glmmTMB/issues
* Date/Publication: 2020-02-03 21:10:05 UTC
* Number of recursive dependencies: 141

Run `revdep_details(,"glmmTMB")` for more info

</details>

## In both

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/AAAtest-all.R’ failed.
    Last 13 lines of output:
      complex bases: ..................
      REML: ......
      Saving and loading glmmTMB objects: .
      ....variance structures: ..............
      weight: .......1.
      ZI models: ..........
      
      ══ Failed ══════════════════════════════════════════════════════════════════════
      ── 1. Failure: Estimates are the same (@test-weight.R#59)  ─────────────────────
      ranef(wei_glmmtmb) not equal to ranef(ind_glmmtmb).
      Component "cond": Component "i": Attributes: < Component "condVar": Mean relative difference: 1.015283e-05 >
      
      ══ DONE ════════════════════════════════════════════════════════════════════════
      Error: Test failures
      Execution halted
    ```

*   checking installed package size ... NOTE
    ```
      installed size is 64.8Mb
      sub-directories of 1Mb or more:
        doc         1.1Mb
        libs       61.3Mb
        test_data   1.2Mb
    ```

# GLMpack

<details>

* Version: 0.1.0
* Source code: https://github.com/cran/GLMpack
* Date/Publication: 2019-07-19 09:40:05 UTC
* Number of recursive dependencies: 81

Run `revdep_details(,"GLMpack")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘AER’ ‘MASS’ ‘Matrix’ ‘censReg’ ‘effects’ ‘foreign’ ‘lme4’ ‘lmtest’
      ‘nnet’ ‘pBrackets’ ‘plm’ ‘pscl’ ‘sandwich’
      All declared Imports should be used.
    ```

# glmulti

<details>

* Version: 1.0.7.1
* Source code: https://github.com/cran/glmulti
* Date/Publication: 2019-04-14 10:45:52 UTC
* Number of recursive dependencies: 16

Run `revdep_details(,"glmulti")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    glmulti,character-character-ANY-ANY: no visible global function
      definition for ‘regsubsets’
    Undefined global functions or variables:
      regsubsets
    ```

# gWQS

<details>

* Version: 2.0.0
* Source code: https://github.com/cran/gWQS
* Date/Publication: 2019-08-27 12:40:02 UTC
* Number of recursive dependencies: 130

Run `revdep_details(,"gWQS")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘dplyr’
      All declared Imports should be used.
    ```

# gWQSRS

<details>

* Version: 1.1.0
* Source code: https://github.com/cran/gWQSRS
* Date/Publication: 2020-01-24 14:30:10 UTC
* Number of recursive dependencies: 104

Run `revdep_details(,"gWQSRS")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘broom’ ‘dplyr’ ‘ggrepel’ ‘kableExtra’ ‘knitr’ ‘nnet’ ‘plotROC’
      All declared Imports should be used.
    ```

# hnp

<details>

* Version: 1.2-6
* Source code: https://github.com/cran/hnp
* Date/Publication: 2018-05-21 15:27:55 UTC
* Number of recursive dependencies: 23

Run `revdep_details(,"hnp")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘glmmADMB’
    ```

# idealstan

<details>

* Version: 0.7.2
* Source code: https://github.com/cran/idealstan
* BugReports: https://github.com/saudiwin/idealstan/issues
* Date/Publication: 2019-07-10 15:00:03 UTC
* Number of recursive dependencies: 109

Run `revdep_details(,"idealstan")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 118.8Mb
      sub-directories of 1Mb or more:
        doc     1.1Mb
        libs  117.0Mb
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# NormalBetaPrime

<details>

* Version: 2.2
* Source code: https://github.com/cran/NormalBetaPrime
* Date/Publication: 2019-01-19 22:40:09 UTC
* Number of recursive dependencies: 13

Run `revdep_details(,"NormalBetaPrime")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  5.1Mb
      sub-directories of 1Mb or more:
        data   4.8Mb
    ```

# parameters

<details>

* Version: 0.5.0
* Source code: https://github.com/cran/parameters
* URL: https://easystats.github.io/parameters/
* BugReports: https://github.com/easystats/parameters/issues
* Date/Publication: 2020-02-09 19:40:03 UTC
* Number of recursive dependencies: 337

Run `revdep_details(,"parameters")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘parameters-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: model_parameters.aov
    > ### Title: Parameters from ANOVAs
    > ### Aliases: model_parameters.aov
    > 
    > ### ** Examples
    > 
    > df <- iris
    > df$Sepal.Big <- ifelse(df$Sepal.Width >= 3, "Yes", "No")
    > 
    > model <- aov(Sepal.Length ~ Sepal.Big, data = df)
    > model_parameters(model, omega_squared = "partial", eta_squared = "partial", epsilon_squared = TRUE)
    Error in model_parameters.aov(model, omega_squared = "partial", eta_squared = "partial",  : 
      Package 'effectsize' required for this function to work. Please install it.
    Calls: model_parameters -> model_parameters.aov
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      ══ testthat results  ═══════════════════════════════════════════════════════════
      [ OK: 285 | SKIPPED: 0 | WARNINGS: 8 | FAILED: 13 ]
      1. Error: efa-cfa (@test-model_parameters.efa_cfa.R#48) 
      2. Failure: model_parameters, standardize-refit (@test-model_parameters_std.R#11) 
      3. Failure: model_parameters, standardize-refit (@test-model_parameters_std.R#12) 
      4. Failure: model_parameters, standardize-refit (@test-model_parameters_std.R#13) 
      5. Failure: model_parameters, standardize-posthoc (@test-model_parameters_std.R#19) 
      6. Failure: model_parameters, standardize-posthoc (@test-model_parameters_std.R#20) 
      7. Failure: model_parameters, standardize-posthoc (@test-model_parameters_std.R#21) 
      8. Failure: model_parameters, standardize-basic (@test-model_parameters_std.R#27) 
      9. Failure: model_parameters, standardize-basic (@test-model_parameters_std.R#28) 
      1. ...
      
      Error: testthat unit tests failed
      Execution halted
    ```

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'effectsize', 'performance'
    ```

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘effectsize’
    ```

# performance

<details>

* Version: 0.4.4
* Source code: https://github.com/cran/performance
* URL: https://easystats.github.io/performance/
* BugReports: https://github.com/easystats/performance/issues
* Date/Publication: 2020-02-10 21:50:06 UTC
* Number of recursive dependencies: 207

Run `revdep_details(,"performance")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘performance-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: binned_residuals
    > ### Title: Binned residuals for logistic regression
    > ### Aliases: binned_residuals
    > 
    > ### ** Examples
    > 
    > model <- glm(vs ~ wt + mpg, data = mtcars, family = "binomial")
    > binned_residuals(model)
    Warning: Probably bad model fit. Only about 50% of the residuals are inside the error bounds.
    Error in print.binned_residuals(x) : 
      Package 'see' required to plot binned residuals. Please install it.
    Calls: <Anonymous> -> print.binned_residuals
    Execution halted
    ```

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘parameters’
    ```

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘parameters’
    ```

# scanstatistics

<details>

* Version: 1.0.1
* Source code: https://github.com/cran/scanstatistics
* URL: https://github.com/BenjaK/scanstatistics
* BugReports: https://github.com/BenjaK/scanstatistics/issues
* Date/Publication: 2018-01-24 12:37:44 UTC
* Number of recursive dependencies: 83

Run `revdep_details(,"scanstatistics")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 11.2Mb
      sub-directories of 1Mb or more:
        libs  10.6Mb
    ```

# sjPlot

<details>

* Version: 2.8.2
* Source code: https://github.com/cran/sjPlot
* URL: https://strengejacke.github.io/sjPlot/
* BugReports: https://github.com/strengejacke/sjPlot/issues
* Date/Publication: 2020-01-23 13:30:02 UTC
* Number of recursive dependencies: 189

Run `revdep_details(,"sjPlot")` for more info

</details>

## In both

*   checking package dependencies ... ERROR
    ```
    Packages required but not available:
      'effectsize', 'ggeffects', 'parameters', 'performance'
    
    See section ‘The DESCRIPTION file’ in the ‘Writing R Extensions’
    manual.
    ```

# sjstats

<details>

* Version: 0.17.9
* Source code: https://github.com/cran/sjstats
* URL: https://github.com/strengejacke/sjstats, https://strengejacke.github.io/sjstats
* BugReports: https://github.com/strengejacke/sjstats/issues
* Date/Publication: 2020-02-06 17:50:02 UTC
* Number of recursive dependencies: 196

Run `revdep_details(,"sjstats")` for more info

</details>

## In both

*   checking package dependencies ... ERROR
    ```
    Packages required but not available:
      'effectsize', 'parameters', 'performance'
    
    Package suggested but not available for checking: ‘sjPlot’
    
    See section ‘The DESCRIPTION file’ in the ‘Writing R Extensions’
    manual.
    ```

