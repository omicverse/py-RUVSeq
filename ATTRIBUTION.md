# Attribution

Upstream RUVSeq 1.33.1, commit f560e84eb493f0c2d72279612f9e4ca99493b434.
Algorithms reconstructed from the upstream source under the accompanying LICENSE.
Upstream author declarations:

```
Package: RUVSeq
Version: 1.33.1
Title: Remove Unwanted Variation from RNA-Seq Data
Description: This package implements the remove unwanted variation
        (RUV) methods of Risso et al. (2014) for the normalization of
        RNA-Seq read counts between samples.
Authors@R: c(person("Davide", "Risso", email = "risso.davide@gmail.com",
	          role = c("aut", "cre", "cph")),
	     person("Sandrine", "Dudoit", role = "aut"),
	     person("Lorena", "Pantano", role = "ctb"),
	     person("Kamil", "Slowikowski", role = "ctb"))
Author: Davide Risso [aut, cre, cph], Sandrine Dudoit [aut], Lorena Pantano [ctb], Kamil Slowikowski [ctb]
Maintainer: Davide Risso <risso.davide@gmail.com>
Date: 04-15-2014
Imports: methods, MASS
Depends: Biobase, EDASeq (>= 1.99.1), edgeR
Suggests: BiocStyle, knitr, RColorBrewer, zebrafishRNASeq, DESeq2
VignetteBuilder: knitr
License: Artistic-2.0
LazyLoad: yes
biocViews: ImmunoOncology, DifferentialExpression, Preprocessing, RNASeq, Software
URL: https://github.com/drisso/RUVSeq
BugReports: https://github.com/drisso/RUVSeq/issues
RoxygenNote: 6.0.1

```
