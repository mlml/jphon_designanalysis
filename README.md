# jphon_designanalysis
Power &amp; effect size error calculations -- materials for Kirby &amp; Sonderegger paper under review


`powerEs.R`: core functions to fit models, run simulations using
`simR` functionality.

`roettgerEtAlData.csv`: dataset from Roettger et al. (2014)
*J. Phonetics*

`runSimsJphon.R`: script to run actual simulations described in article, resulting in RDS
file. (Warning: takes several days!)

`summarizeVisualizeRuns.R`: script to generate plots from submitted
article.

`pwr examples.Rmd`: examples of using `pwr` package for simple power
calculations mentioned in article text.


To actually run simulations generated in paper:

* `source(runSimsJphon.R)`.

Because this takes several days, the resulting data file is inlcluded
(`runs_26Sep17_nRuns500.rds`), and you can do data summary and
visualizaton reported in the paper by just running:

* `source(summarizeVisualizeRuns.R`
