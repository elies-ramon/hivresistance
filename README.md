# hivresistance

This R code allows to reproduce the experiments in the paper: 

Ramon, Elies. *Unraveling HIV protease drug resistance and genetic diversity with kernel methods.* bioRxiv 2025.03.26.644092; doi: [https://doi.org/10.1101/2025.03.26.644092](https://www.biorxiv.org/content/10.1101/2025.03.26.644092v1).

The required packages are:

* [kerntools](https://cran.r-project.org/web/packages/kerntools/index.html)
* kernlab
* maotai
* stringi

(Plots)
* ggplot2
* reshape2
* cowplot
* viridis


The workflow is:

1. preprocessing.R
2. residue_dist.R
3. kernel_matrices.R
4. splitting.R
5. linkage_diseq.R
6. descriptive.R
7. kpca.R
8. mixtures_vs_nomixtures.R
9. svm.R
10. results.R

If you detect any error or want to ask anything about this work, you can open an issue or e-mail me at: eramon@everlyrusher.com.
