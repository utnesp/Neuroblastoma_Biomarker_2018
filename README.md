# Neuroblastoma_Biomarker_2018
Scripts used in "Clinically relevant biomarker discovery in high-risk recurrent neuroblastoma"

SEQC.SampleCharacteristics.txt was generated from the series matrix files associated with Gene Expression Omnibuss accession [GSE62564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62564) and [GSE49711](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711) using the R script [getSeriesMatrixCharacteristics.R](https://github.com/utnesp/NORAD/blob/master/getSeriesMatrixCharacteristics.R)

Scatterplot regarding RNA-seq validaton was generated based on qPCR data from the two files:
- RTqPCR_validation.signal.txt
- RTqPCR_validation.signal.sdev.txt
The data was plotted using ggplot2 with the following command:
```r
ggplot(data, aes(qPCR, SEQ)) + geom_point() + facet_wrap(~external_gene_name, scales = "free", labeller=labeller(external_gene_name = unlist(gene_names))) + theme_bw()
```

PCA and PLS-DA result were generated using PCA_PLSDA.R 

NBCellLinesRawCounts.txt represent raw counts in neuroblastoma cell lines used in this study.

## Characterization of individual RNA classes
An interactive plot characterizing individual RNA classes is rendered as a html through [GitHub pages](https://utnesp.github.io/Neuroblastoma_Biomarker_2018/)

Tutorials on how to create websites with R were used as guidelines, including:
[R Markdown](https://rmarkdown.rstudio.com/rmarkdown_websites.html)
[GitHub pages](https://pages.github.com)

The following files were used to deploy the rendered html: 
- expressionRNAclasses.Rmd (used to generated index.html)
- [docs/CPM.cumulative.freq.txt](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/docs/CPM.cumulative.freq.txt) (used to generated interactive plot within expressionRNAclasses.Rmd)
- index.html (knitted within RStudio from expressionRNAclasses.Rmd)
- index.Rmd (defines shared output options for all R Markdown documents within a site)
- _config.yml (configuration file)
- _site.yml (defines common navigation bar)



