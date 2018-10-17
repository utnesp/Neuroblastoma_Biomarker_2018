# Neuroblastoma_Biomarker_2018
Scripts used in "Clinically relevant biomarker discovery in high-risk recurrent neuroblastoma"

- SEQC.SampleCharacteristics.txt was generated from the series matrix files associated with Gene Expression Omnibuss accession [GSE62564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62564) and [GSE49711](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711) using the R script [getSeriesMatrixCharacteristics.R](https://github.com/utnesp/NORAD/blob/master/getSeriesMatrixCharacteristics.R)

Figure regarding RNA-seq validaton was generated using ggplot2 from raw data files:
- RTqPCR_validation.signal.txt
- RTqPCR_validation.signal.sdev.txt

PCA and PLS-DA result was geneated using PCA_PLSDA.R 

NBCellLinesRawCounts.txt represent raw counts in neuroblastoma cell lines used in this study.

# Characterization of individual RNA classes
An interactive plot characterizing individual RNA classes is rendered as a html through [GitHub pages](https://utnesp.github.io/Neuroblastoma_Biomarker_2018/)

Tutorials on how to create websites with R were used as guidelines, and include:
[R Markdown](https://rmarkdown.rstudio.com/rmarkdown_websites.html)
[GitHub pages](https://pages.github.com)

The following files were used to deploy the rendered html. 
expressionRNAclasses.Rmd (used to generated index.html)
index.html (knitted from expressionRNAclasses.Rmd)
index.Rmd (defines shared output options for all R Markdown documents within a site)
_config.yml (configuration file)
_site.yml (defines common navigation bar)



