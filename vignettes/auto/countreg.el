(TeX-add-style-hook "countreg"
 (lambda ()
    (LaTeX-add-bibliographies)
    (LaTeX-add-labels
     "sec:intro"
     "sec:software"
     "tab:overview"
     "eq:family"
     "eq:mean"
     "eq:Poisson"
     "eq:negbin"
     "eq:hurdle"
     "eq:hurdle-mean"
     "eq:zeroinfl"
     "eq:zeroinfl-mean"
     "sec:illustrations"
     "fig:ofp"
     "fig:bad-good"
     "fig:ofp2"
     "tab:summary"
     "sec:summary"
     "app:hurdle"
     "app:zeroinfl"
     "app:methods"
     "tab:methods"
     "app:replication")
    (TeX-add-symbols
     '("fct" 1)
     '("class" 1))
    (TeX-run-style-hooks
     "thumbpdf"
     ""
     "latex2e"
     "jss10"
     "jss"
     "nojss")))

