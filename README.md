# Neuroblastoma Biomarker 2018
Scripts used in "Clinically relevant biomarker discovery in high-risk recurrent neuroblastoma"

## Generation of counts
Fastq files sequenced on the Illumina HiSeq 2000 system were deposited to the Sequence Read Archive with accession [PRJNA491629](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA491629). Fastq files were further mapped to the human genome build GRCh38 using [STAR aligner](https://github.com/alexdobin/STAR) and [STAR.sh](https://github.com/utnesp/Neuroblastoma_Biomarker_2018/blob/master/STAR.sh). [featureCounts](http://subread.sourceforge.net/) and the script [featureCounts.sh](https://github.com/utnesp/Neuroblastoma_Biomarker_2018/blob/master/featureCounts.sh) was used to generate [NBCellLinesRawCounts.txt](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/NBCellLinesRawCounts.txt) which represent raw counts in neuroblastoma cell lines used in this study. The human genome sequence (primary assembly) was obtained from [GENCODE](https://www.gencodegenes.org/human/), whereas genome annotation was obtained from [ENSEMBL](http://www.ensembl.org/info/data/ftp/index.html).

## Generation of SEQC sample characteristics
[SEQC.SampleCharacteristics.txt](https://github.com/utnesp/Neuroblastoma_Biomarker_2018/blob/master/SEQC.SampleCharacteristics.txt) was generated from the series matrix files associated with Gene Expression Omnibuss accession [GSE62564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62564) and [GSE49711](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711) using the R script [getSeriesMatrixCharacteristics.R](https://github.com/utnesp/NORAD/blob/master/getSeriesMatrixCharacteristics.R) and the R snippet:

```r
ftpURLof_series_matrix.txt.gz = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62564/matrix/GSE62564_series_matrix.txt.gz" # reanalyzed dataset
series_matrix <- getSeriesMatrixCharacteristics(ftpURLof_series_matrix.txt.gz)
ftpURLof_series_matrix.txt.gz = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49711/matrix/GSE49711_series_matrix.txt.gz" # original dataset (also includes age_at_diagnosis, mycn_status, inss_stage and class_label)
series_matrix2 <- getSeriesMatrixCharacteristics(ftpURLof_series_matrix.txt.gz)
SEQC.ClincicalAttributes <- cbind(series_matrix, series_matrix2[, c("age_at_diagnosis", "mycn_status", "inss_stage", "class_label")]) # see descriptions at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711
```

The resulting file may be utilized within R using the data.table command fread:
```r
SEQC.ClincicalAttributes <- fread("https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/SEQC.SampleCharacteristics.txt")
```

## RNA-seq validation
Scatterplot regarding RNA-seq validaton was generated based on normalized qPCR data from the two files:
- [RTqPCR_validation.signal.txt](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/RTqPCR_validation.signal.txt)
- [RTqPCR_validation.signal.sdev.txt](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/RTqPCR_validation.signal.sdev.txt)

The two files were melted to long format (see [Introduction to reshape2](https://seananderson.ca/2013/10/19/reshape/)) and the resulting data with the two new columns qPCR and SEQ was plotted using [ggplot2](https://github.com/tidyverse/ggplot2) with the following command:
```r
ggplot(data, aes(qPCR, SEQ)) + geom_point() + facet_wrap(~external_gene_name, scales = "free", labeller=labeller(external_gene_name = unlist(gene_names))) + theme_bw()
```
Correlation tests were carried out using the R base cor.test() function.

## PCA PLS-DA
PCA and PLS-DA results were generated using [PCA_PLSDA.R](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/PCA_PLSDA.R) following [tutorials using the mixOmics package](http://mixomics.org/case-studies/).


## Characterization of individual RNA classes
An interactive plot characterizing individual RNA classes is rendered as a html through [GitHub pages](https://utnesp.github.io/Neuroblastoma_Biomarker_2018/).

Tutorials on how to create websites with R were used as guidelines, including:
- [R Markdown](https://rmarkdown.rstudio.com/rmarkdown_websites.html)
- [GitHub pages](https://pages.github.com)

The following files were used to deploy the rendered html: 
- [expressionRNAclasses.Rmd](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/expressionRNAclasses.Rmd)(used to generated index.html)
- [docs/CPM.cumulative.freq.txt](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/docs/CPM.cumulative.freq.txt) (used to generate interactive plot within expressionRNAclasses.Rmd)
- [index.html](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/index.html) (knitted within RStudio from expressionRNAclasses.Rmd)
- [index.Rmd](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/index.Rmd) (defines shared output options for all R Markdown documents within a site)
- [_config.yml](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/_config.yml) (configuration file)
- [_site.yml](https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/_site.yml) (defines common navigation bar)


Alternatively, running this code snippet will also generate an interactive plot:

```r
library(ggplot2)
library(plotly)
library(shiny)
library(data.table)

data <- fread("https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/docs/CPM.cumulative.freq.txt")

shapes = c(15, 16, 17, 0, 18, 3, 5, 2, 8, 24, 1, 25, 4, 7, 8, 9, 10, 26:30)
s <- data.frame(gene_biotype = data$gene_biotype[!duplicated(data$gene_biotype)], shape = shapes[1:length(unique(data$col))] , color = unique(data$col))
p1 <- suppressWarnings(ggplot(data, aes(CPM, CumulativeFrequency)) + geom_point(aes(shape = gene_biotype, color = gene_biotype, 
 text = paste("GENE: ", paste(external_gene_name, "(", ensembl_gene_id, ")", sep = ""), "\nGENE LENGTH (bp):", GeneLength)), size = 1) + theme_classic() + labs(x = "Average expression across cell lines (log2 CPM)", y = "Cumulative Relative Frequency", color='Biotype') + theme(legend.position = "right") + 
       scale_colour_manual(name = "Biotype",
                            labels = s$gene_biotype,
                            values = as.character(s$color))  +
        scale_shape_manual(name = "Biotype",
                           labels = s$gene_biotype,
                           values = s$shape) +
        theme(text = element_text(size = 10), axis.text = element_text(size = 10)))
ggplotly(p1 + theme(text = element_text(size = 10), axis.text = element_text(size = 10), legend.text = element_text(size = 10)), width = 1000, height = 1000)
