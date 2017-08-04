# Setup

## Platform

|setting  |value                        |
|:--------|:----------------------------|
|version  |R version 3.4.1 (2017-06-30) |
|system   |x86_64, darwin15.6.0         |
|ui       |RStudio (1.0.153)            |
|language |(EN)                         |
|collate  |en_AU.UTF-8                  |
|tz       |Australia/Sydney             |
|date     |2017-08-03                   |

## Packages

|package |*  |version |date       |source                |
|:-------|:--|:-------|:----------|:---------------------|
|pscl    |   |1.5.1   |2017-08-03 |local (atahk/pscl@NA) |

# Check results

30 packages

|package        |version | errors| warnings| notes|
|:--------------|:-------|------:|--------:|-----:|
|accelmissing   |1.1     |      0|        0|     0|
|AER            |1.2-5   |      0|        0|     0|
|agridat        |1.12    |      0|        0|     0|
|AICcmodavg     |2.1-1   |      0|        0|     0|
|anominate      |0.5     |      1|        0|     0|
|ashr           |2.0.5   |      0|        0|     2|
|BANFF          |2.0     |      0|        0|     0|
|BayesSpec      |0.5.3   |      0|        0|     0|
|BSGS           |2.0     |      0|        0|     1|
|catdata        |1.2.1   |      0|        0|     0|
|congressbr     |0.1.1   |      0|        0|     0|
|DAMisc         |1.4-3   |      0|        0|     0|
|DClusterm      |0.1     |      0|        0|     1|
|emIRT          |0.0.8   |      1|        0|     0|
|fakeR          |1.0     |      0|        0|     0|
|fence          |1.0     |      0|        0|     0|
|glmmTMB        |0.1.1   |      0|        0|     1|
|glmulti        |1.0.7   |      1|        0|     0|
|glmx           |0.1-1   |      0|        0|     0|
|hnp            |1.2-3   |      0|        0|     1|
|lctools        |0.2-5   |      0|        0|     0|
|lsmeans        |2.26-3  |      0|        0|     1|
|mpath          |0.2-4   |      1|        0|     0|
|oc             |0.96    |      1|        0|     0|
|opticut        |0.1-0   |      0|        0|     0|
|paleofire      |1.1.9   |      0|        0|     0|
|robustsae      |0.1.0   |      0|        0|     0|
|sandwich       |2.4-0   |      0|        0|     0|
|scanstatistics |0.1.0   |      0|        0|     0|
|wnominate      |1.2     |      1|        0|     0|

## accelmissing (1.1)
Maintainer: Jung Ae Lee <jungaeleeb@gmail.com>

0 errors | 0 warnings | 0 notes

## AER (1.2-5)
Maintainer: Achim Zeileis <Achim.Zeileis@R-project.org>

0 errors | 0 warnings | 0 notes

## agridat (1.12)
Maintainer: Kevin Wright <kw.stat@gmail.com>  
Bug reports: https://github.com/kwstat/agridat/issues

0 errors | 0 warnings | 0 notes

## AICcmodavg (2.1-1)
Maintainer: Marc J. Mazerolle <marc.mazerolle@sbf.ulaval.ca>

0 errors | 0 warnings | 0 notes

## anominate (0.5)
Maintainer: Christopher Hare <chare@uga.edu>

1 error  | 0 warnings | 0 notes

```
checking whether package ‘anominate’ can be installed ... ERROR
Installation failed.
See ‘/Users/zoemeers/pscl_2/revdep/checks/anominate.Rcheck/00install.out’ for details.
```

## ashr (2.0.5)
Maintainer: Peter Carbonetto <pcarbo@uchicago.edu>

0 errors | 0 warnings | 2 notes

```
checking package dependencies ... NOTE
Packages which this enhances but not available for checking:
  ‘REBayes’ ‘Rmosek’

checking compiled code ... NOTE
File ‘ashr/libs/ashr.so’:
  Found no calls to: ‘R_registerRoutines’, ‘R_useDynamicSymbols’

It is good practice to register native routines and to disable symbol
search.

See ‘Writing portable packages’ in the ‘Writing R Extensions’ manual.
```

## BANFF (2.0)
Maintainer: Tianwei Yu <tianwei.yu@emory.edu>

0 errors | 0 warnings | 0 notes

## BayesSpec (0.5.3)
Maintainer: Andrew Ferris <andrew.ferris@sydney.edu.au>

0 errors | 0 warnings | 0 notes

## BSGS (2.0)
Maintainer: Kuo-Jung Lee <kuojunglee@mail.ncku.edu.tw>

0 errors | 0 warnings | 1 note 

```
checking R code for possible problems ... NOTE
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

## catdata (1.2.1)
Maintainer: Gunther Schauberger
 <gunther.schauberger@stat.uni-muenchen.de>

0 errors | 0 warnings | 0 notes

## congressbr (0.1.1)
Maintainer: Robert Myles McDonnell <robertmylesmcdonnell@gmail.com>

0 errors | 0 warnings | 0 notes

## DAMisc (1.4-3)
Maintainer: Dave Armstrong <dave@quantoid.net>

0 errors | 0 warnings | 0 notes

## DClusterm (0.1)
Maintainer: Virgilio Gomez-Rubio <virgilio.gomez@uclm.es>

0 errors | 0 warnings | 1 note 

```
checking package dependencies ... NOTE
Package suggested but not available for checking: ‘INLA’
```

## emIRT (0.0.8)
Maintainer: James Lo <lojames@usc.edu>

1 error  | 0 warnings | 0 notes

```
checking whether package ‘emIRT’ can be installed ... ERROR
Installation failed.
See ‘/Users/zoemeers/pscl_2/revdep/checks/emIRT.Rcheck/00install.out’ for details.
```

## fakeR (1.0)
Maintainer: Lily Zhang <lilyhzhang1029@gmail.com>

0 errors | 0 warnings | 0 notes

## fence (1.0)
Maintainer: Thuan Nguyen <nguythua@ohsu.edu>

0 errors | 0 warnings | 0 notes

## glmmTMB (0.1.1)
Maintainer: Mollie Brooks <mollieebrooks@gmail.com>  
Bug reports: https://github.com/glmmTMB/glmmTMB/issues

0 errors | 0 warnings | 1 note 

```
checking installed package size ... NOTE
  installed size is 15.2Mb
  sub-directories of 1Mb or more:
    doc    1.2Mb
    libs  13.5Mb
```

## glmulti (1.0.7)
Maintainer: Vincent Calcagno <vincent.calcagno@sophia.inra.fr>

1 error  | 0 warnings | 0 notes

```
checking whether package ‘glmulti’ can be installed ... ERROR
Installation failed.
See ‘/Users/zoemeers/pscl_2/revdep/checks/glmulti.Rcheck/00install.out’ for details.
```

## glmx (0.1-1)
Maintainer: Achim Zeileis <Achim.Zeileis@R-project.org>

0 errors | 0 warnings | 0 notes

## hnp (1.2-3)
Maintainer: Rafael de Andrade Moral <rafael_moral@yahoo.com.br>

0 errors | 0 warnings | 1 note 

```
checking package dependencies ... NOTE
Package suggested but not available for checking: ‘glmmADMB’
```

## lctools (0.2-5)
Maintainer: Stamatis Kalogirou <skalo@hua.gr>

0 errors | 0 warnings | 0 notes

## lsmeans (2.26-3)
Maintainer: Russell Lenth <russell-lenth@uiowa.edu>  
Bug reports: https://github.com/rvlenth/lsmeans/issues

0 errors | 0 warnings | 1 note 

```
checking package dependencies ... NOTE
Packages suggested but not available for checking: ‘glmmADMB’ ‘lme4.0’
```

## mpath (0.2-4)
Maintainer: Zhu Wang <zwang@connecticutchildrens.org>

1 error  | 0 warnings | 0 notes

```
checking whether package ‘mpath’ can be installed ... ERROR
Installation failed.
See ‘/Users/zoemeers/pscl_2/revdep/checks/mpath.Rcheck/00install.out’ for details.
```

## oc (0.96)
Maintainer: James Lo <lojames@usc.edu>

1 error  | 0 warnings | 0 notes

```
checking whether package ‘oc’ can be installed ... ERROR
Installation failed.
See ‘/Users/zoemeers/pscl_2/revdep/checks/oc.Rcheck/00install.out’ for details.
```

## opticut (0.1-0)
Maintainer: Peter Solymos <solymos@ualberta.ca>  
Bug reports: https://github.com/psolymos/opticut/issues

0 errors | 0 warnings | 0 notes

## paleofire (1.1.9)
Maintainer: Olivier Blarquez <blarquez@gmail.com>

0 errors | 0 warnings | 0 notes

## robustsae (0.1.0)
Maintainer: Jiyoun Myung <jiyoun@ufl.edu>

0 errors | 0 warnings | 0 notes

## sandwich (2.4-0)
Maintainer: Achim Zeileis <Achim.Zeileis@R-project.org>

0 errors | 0 warnings | 0 notes

## scanstatistics (0.1.0)
Maintainer: Benjamin Kjellson <benjak@math.su.se>  
Bug reports: https://github.com/BenjaK/scanstatistics/issues

0 errors | 0 warnings | 0 notes

## wnominate (1.2)
Maintainer: James Lo <lojames@usc.edu>

1 error  | 0 warnings | 0 notes

```
checking whether package ‘wnominate’ can be installed ... ERROR
Installation failed.
See ‘/Users/zoemeers/pscl_2/revdep/checks/wnominate.Rcheck/00install.out’ for details.
```

