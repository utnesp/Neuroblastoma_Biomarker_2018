# load required packages
library(Biobase)
library(biomaRt)
library(clusterProfiler)
library(cum.plot) # To install: devtools::install_github("utnesp/cumulative.plot.R")
library(data.table)
library(easybiomart) # To install: # devtools::install_github("utnesp/Easy-bioMart")
library(forestplot)
library(gplots)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(grid)
library(gridExtra)
library(gridBase)
library(gridGraphics)
library(lattice)
library(metagenomeSeq)
library(mixOmics)
library(org.Hs.eg.db)
library(pheatmap)
library(plot3D)
library(plotly)
library(RColorBrewer)
library(RCurl)
library(ReactomePA)
library(reshape2)
library(Rgraphviz)
library(RISmed)
library(S4Vectors) 
library(scatterplot3d)
library(splitstackshape)
library(stringr) 
library(survival)
library(survminer)

source_https <- function(url, ...) {
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

source_https("https://raw.githubusercontent.com/utnesp/functions/master/R/functions.R")

setEPS() # to be able to save as eps
# set annotation database
mart = useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl')

# set working directory
folder = "~/Desktop/test" # change the name of the folder to your choosing
dir.create(folder,recursive = T) # create dir if needed
setwd(folder)  # set working directory 

# read in raw expression
PrePostCellLines <- fread("https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/NBCellLinesRawCounts.txt")

# create design with repeated measures
design <- data.frame(body_site_s = c(rep("Primary Tumor", 5), rep("Recurrent Tumor", 5)), id_s = rep(c("BE2", "KANR", "NBLWR", "CHLA136", "CHLA20"), 2))

# normalize samples based on protein coding genes before using mixOmics
biotypes <- ensg2ext_name_biotype(PrePostCellLines$ensembl_gene_id)
y <- DGEList(PrePostCellLines[, 2:ncol(PrePostCellLines)], genes = PrePostCellLines$ensembl_gene_id)
temp <- PrePostCellLines[PrePostCellLines$ensembl_gene_id %in% biotypes[biotypes$gene_biotype == "protein_coding", "ensembl_gene_id"], ]
y.norm <- DGEList(temp[, 2:ncol(temp)], genes = temp$ensembl_gene_id) # normalize only on protein coding content
y.norm <- edgeR::calcNormFactors(y.norm, method = "RLE")
y$samples <- y.norm$samples

X <- cpm(y)
row.names(X) <- y$genes$genes

### check if RNA-seq validates with qPCR 
validation <- fread("https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/RTqPCR_validation.signal.txt")
validation.sd <- fread("https://raw.githubusercontent.com/utnesp/Neuroblastoma_Biomarker_2018/master/RTqPCR_validation.signal.sdev.txt")
counts <- X
data <- melt(validation, id.vars = "gene_name")
data.sd <- melt(validation.sd, id.vars = "gene_name")
data$sd <- data.sd$value
primers = data.frame(external_gene_name = c("FN1", "GAS5", "HDAC1", "MALAT1", "MYCN", "MYCNOS", "NEAT1", "SNHG16", "TFRC"))
primers <- hugo2ensg(primers$external_gene_name, combine = T)
counts <- counts[row.names(counts) %in% primers$ensembl_gene_id, ]
names <- ensg2ext_name(row.names(counts))
row.names(counts) <- names[match(row.names(counts), names$ensembl_gene_id) ,"external_gene_name"]
counts <- as.data.frame(counts)
counts <- log2(counts+1)
counts$external_gene_name <- row.names(counts)
data2 <- melt(counts, id.vars = "external_gene_name")
colnames(data2)[3] <- "SEQ" 
colnames(data)[3] <- "qPCR"
colnames(data)[1] <- colnames(data2)[1]
data <- merge(data,data2)
data$external_gene_name <- as.factor(data$external_gene_name)
cordata <- data.frame(gene = levels(data$external_gene_name))
cordata$cor <- NA
    for (i in 1:length(levels(data$external_gene_name))) {
        subdata <- data[data$external_gene_name == levels(data$external_gene_name)[i], ]
        summary <- cor.test(subdata$qPCR, subdata$SEQ)
        cordata$cor[i] <- round(summary$estimate,2)
        cordata$pvalue[i] <- signif(summary$p.value,2)
    }

cordata$gene <- factor(cordata$gene, levels = levels(data$external_gene_name))
data <- data[data$external_gene_name %in% cordata$gene, ]
gene_names <- as.list(paste(cordata$gene, "\nr=", cordata$cor, " p=", cordata$pvalue, sep = ""))
names(gene_names) <- cordata$gene
#postscript("Figure1.eps") # uncomment to save
ggplot(data, aes(qPCR, SEQ)) + geom_point() + facet_wrap(~external_gene_name, scales = "free", labeller=labeller(external_gene_name = unlist(gene_names))) + theme_bw() + labs(x= "Relative expression compared to reference genes (RT-qPCR)", y="log2 CPM (RNA-seq)")
#dev.off() # uncomment to save

X <- X+1 # add priour count of 1 to avoid log2 errors
# convert X to integers for faster computations (for mixOmics analysis)
X <- apply(X, c(1,2), function(X) as.integer(round(X, 0)))
X <- t(X) # transpose matrix (required by mixOmics)

# use gene names as row.names instead of gene ids (for visualization purposes)
biotypes <- ensg2ext_name_biotype(colnames(X))
dups <- biotypes[duplicated(biotypes$external_gene_name), ] # genes that have duplicated names / corresponding to more than one gene id 
dups$external_gene_name <- paste(dups$external_gene_name, dups$ensembl_gene_id, sep = ".") # append gene id to them to make them unique
biotypes$external_gene_name[which(duplicated(biotypes$external_gene_name))] <- dups$external_gene_name
X <- X[, colnames(X) %in% biotypes$ensembl_gene_id] 
biotypes <- biotypes[match(colnames(X), biotypes$ensembl_gene_id), ]
identical(biotypes$ensembl_gene_id, colnames(X))
x <- X # keep a list with ensg identifiers for later analysis on biotypes
colnames(X) <- biotypes$external_gene_name 
row.names(design) <- row.names(X)

# Keep only genes that have counts more than threshold in at least 4 samples in each group
threshold = 10
Z <- apply(X, c(1,2), function(x) x > threshold)
# keep only genes in primary tumors where at least 4 samples have count more than threshold
G <- apply(Z[1:5, ], 2, function(x) sum(as.integer(x)))
idx.pri <- which(as.numeric(G) >= 4)
# keep only genes in recurrent tumors where at least 4 samples have count more than threshold
G <- apply(Z[6:10, ], 2, function(x) sum(as.integer(x)))
idx.rec <- which(as.numeric(G) >= 4)
# display how many were unique to primary and to recurrent tumors, and how many were shared
input <- list(primary = idx.pri, recurrent = idx.rec)
idx <- c(idx.pri, idx.rec)
idx <- idx[order(idx)]
idx <- idx[!duplicated(idx)]
X <- X[, idx]
x <- x[, idx]

# split genes on proteins and lncRNA
lnc <- biotypes[biotypes$gene_biotype != "protein_coding", "external_gene_name"]
X.lnc <- X[, colnames(X) %in% lnc]
X <- X[, !colnames(X) %in% lnc] 

# PCA
pca <- pca(X, ncomp = 5, logratio = 'CLR', multilevel = design$id_s, scale =F)
# PCA
p <- data.table(`Explained Variance` = pca$explained_variance, `Principal Components` = 1:length(pca$explained_variance))
p <- ggplot(p) + geom_col(aes_(x = p$`Principal Components`, as.numeric(p$`Explained Variance`))) + theme_light() + labs(x = "Principal Components", y = "Explained Variance") + scale_x_continuous(breaks = 1:10)
print(p)

# all cell lines except NBLW separates in PC1
plotIndiv(pca, ind.names = row.names(X),
          group = design$body_site_s,
          title = 'PC1 vs PC2', ellipse = F, legend = F)

# 3d pca
pca.variates <- as.data.frame(pca$variates$X[, 1:3])
# create distance effect in PCA plot (size of dots increase or decrease with shorter or longer distance, resepectively)
coord = apply(pca.variates, 2, min) 
coord <- apply(pca.variates, 1, function(x) x - coord) 
coord = round(apply(coord, 2, function(x) sum(x^2)), 0)
coord = signif(coord, 2)
coord = coord / max(coord)
coord = 1 - coord
coord = ifelse(coord < 0.1, 0.1, coord)
angle = 130; scale.y = 0.8 # which perspective to look at 3d box

s3d.plot <- function() {
s3d <- scatterplot3d(pca.variates$PC1, pca.variates$PC2, pca.variates$PC3, color = color.mixo(design$body_site_s), angle = angle,box = F, cex.symbols = coord * 8,
       xlab = paste('PC1:', round(pca$explained_variance,2)[1]*100, '% expl. var'),
       ylab = paste('PC2:', round(pca$explained_variance,2)[2]*100, '% expl. var'),
       zlab = paste('PC3:', round(pca$explained_variance,2)[3]*100, '% expl. var'), 
       label.tick.marks = T, scale.y = scale.y, highlight.3d = F,
        pch = 19,
       surface=FALSE, grid = T,
       colkey = F, contour3d = T,
       axis.scales = FALSE, labels = row.names(pca.variates), id.n=nrow(pca.variates)
       )
txt.pos = pca.variates[, 1:3]
text(s3d$xyz.convert(txt.pos), labels = rownames(pca.variates),
     cex= 1, col = "black", pos = 4)
title("PCA")
box("figure", lty = 1)
}

grab_grob <- function(){
  grid.echo()
  grid.grab()
}

s3d.plot()
a <- grab_grob() # for use with pls-da plotting 

# 3d biplot 
{
    pca.variates <- pca$variates$X[, 1:3]
    pca.loadings <- pca$loadings$X[,1:3]
    pca.loadings <- as.data.frame(apply(pca.loadings, 2, function(x) x/max(abs(x))))
    pca.variates <- as.data.frame(apply(pca.variates, 2, function(x) x/max(abs(x))))
    pca.variates$names <- row.names(pca.variates)
    pca.variates$color <- color.mixo(design$body_site_s)
    pca.variates$size <- 20
    pca.variates$type = design$body_site_s
    pca.loadings$names <- ""
    pca.loadings$color <- "#999999"
    pca.loadings$size <- (sqrt(pca.loadings$PC1^2 + pca.loadings$PC2^2 + pca.loadings$PC3^2)*10)+1
    pca.loadings$size <- round(pca.loadings$size)
    pca.loadings$type = factor("gene")
    pca.variates <- rbind(pca.variates, pca.loadings)
    pca.variates$rownames <- row.names(pca.variates)
    pca.variates$color <- factor(pca.variates$color, levels = c("#388ECC","#F68B33","#999999"))
    
    p <- plot_ly(pca.variates, x = ~PC1, y = ~PC2, z = ~PC3, mode = 'markers', type = 'scatter3d', color = ~type, size = ~size, colors = levels(pca.variates$color), marker = list(symbol = 'circle', sizemode = 'diameter', size = ~size, text = ~names), sizes = c(min(pca.loadings$size), max(pca.loadings$size)), text=~rownames, hoverinfo  = 'text') %>%
      layout(scene = list(xaxis = list(title = paste('PC1:', round(pca$explained_variance,2)[1]*100, '% expl. var')),
                         yaxis = list(title = paste('PC2:', round(pca$explained_variance,2)[2]*100, '% expl. var')),
                         zaxis = list(title = paste('PC3:', round(pca$explained_variance,2)[3]*100, '% expl. var'))))
    p
}

# PLS-DA
plsda <- plsda(X, Y = design$body_site_s, multilevel = design$id_s, ncomp = 3, logratio = 'CLR', scale = F)
textsize=12
# plot PCA and PLS-DA side-by-side
p1 <- plotIndiv(plsda, ind.names = row.names(X), # data for plotting in plsda$variates$X
          group = design$body_site_s, comp = c(1,2),
          title = 'PLS-DA', ellipse = F, legend = T,  
          , cex = 5, size.title = textsize, size.axis = textsize, size.legend = textsize, size.xlabel = textsize, size.ylabel = textsize, legend.title = "", legend.position = "bottom")
p2 <- as.data.frame(p1$graph$data)
p2$col <- factor(p2$col, levels = c("#F68B33", "#388ECC"))
p2$group <- factor(p2$group, levels = c("Primary Tumor", "Recurrent Tumor"))
cols <- c("Primary Tumor" = "#388ECC", "Recurrent Tumor" = "#F68B33")
p3 <- ggplot(p2, aes(x,y, color = group)) + theme_bw() + geom_point(size = 4) + geom_text_repel(aes(label = names), color = "black", size = 5.5, box.padding = unit(2, "mm")) + scale_color_manual(values = cols, name = "") + ggtitle("PLS-DA") + labs(x = p1$graph$labels$x, y = p1$graph$labels$y) + theme(text = element_text(size = 12), legend.position = "top", plot.title = element_text(hjust = 0.5)) 
grid.arrange(a,p3, ncol = 2) 

# sPLS-DA
ncomp = 2; keepX <- c(300,300)
splsda <- splsda(X, Y = design$body_site_s,multilevel = design$id_s,ncomp = ncomp, keepX = keepX, logratio = 'CLR', scale = F)

# verify that separation in comp1 is still present
plotIndiv(splsda, ind.names = row.names(X), # data for plotting in plsda$variates$X
          group = design$body_site_s,
          title = 'multilevel sPLS-DA, comp 1 - 2', ellipse = F, legend = T,
          , cex = 10, size.title = 20, size.axis = 20, size.legend = 20, size.xlabel = 20, size.ylabel = 20, legend.title = "", legend.position = "bottom")

# cross-validation
perf.splsda <- perf(splsda, validation = "loo", multilevel = design$id_s,
                   logratio = 'CLR', scale = F, progressBar = T, cpus = 7) 
# match the selected variables to the stable features
ind.match = match(selectVar(splsda, comp = 1)$name, 
                  names(perf.splsda$features$stable[[1]]))
#extract the frequency of selection of those selected variables
Freq = as.numeric(perf.splsda$features$stable[[1]][ind.match])
selectedVars <- data.frame(loading = selectVar(splsda, comp = 1)$value, stability = Freq)
stability.cell <- selectedVars[selectedVars$stability >= 1 & !is.na(selectedVars$stability), ]
colnames(stability.cell)[1] <- "loadings.cell.comp1"
colnames(stability.cell)[2] <- "stability.cell.comp1"
stability.cell$external_gene_name <- row.names(stability.cell) 
stability.cell <- ext_name2ensg(stability.cell$external_gene_name, combine = T)
stability.cell <- stability.cell[order(-abs(stability.cell$loadings.cell.comp1)), ]

# plot genes separating the two groups
cutoff = 0.02
pheat <- X[, colnames(X) %in% stability.cell[abs(stability.cell$loadings.cell.comp1) >= cutoff, "external_gene_name"]]
# Keep only genes that have coordinated expression in majority of paired samples (3/5)
log.fc.threshold = 0.5
range.threshold = 10
nr.samples = 3
Z <- apply(log2(pheat[6:10, ] / pheat[1:5, ]), c(1,2), function(x) x >= log.fc.threshold)
# up reg
G <- apply(Z, 2, function(x) sum(as.integer(x)))
idx.pri <- which(as.numeric(G) >= nr.samples)
Z <- apply(pheat[6:10, idx.pri] - pheat[1:5,  idx.pri], c(1,2), function(x) x >= range.threshold)
G <- apply(Z, 2, function(x) sum(as.integer(x)))
idx.pri <- idx.pri[which(as.numeric(G) >= nr.samples)]
# down reg
Z <- apply(log2(pheat[6:10, ] / pheat[1:5, ]), c(1,2), function(x) x <= -log.fc.threshold)
G <- apply(Z, 2, function(x) sum(as.integer(x)))
idx.rec <- which(as.numeric(G) >= nr.samples)
Z <- apply(pheat[6:10, idx.rec] - pheat[1:5,  idx.rec], c(1,2), function(x) x <= -range.threshold)
G <- apply(Z, 2, function(x) sum(as.integer(x)))
idx.rec <- idx.rec[which(as.numeric(G) >= nr.samples)]
# display how many were unique to primary and to recurrent tumors
input <- list(primary = idx.pri, recurrent = idx.rec)
idx <- c(idx.pri, idx.rec)
idx <- idx[order(idx)]
idx <- idx[!duplicated(idx)]
pheat <- pheat[, idx] 
stability.cell$loadings.cell.comp1 <- round(stability.cell$loadings.cell.comp1, 2)
stability.cell <- stability.cell[order(stability.cell$loadings.cell.comp1), ]
ann <- data.frame(Loadings = stability.cell[match(colnames(pheat), stability.cell$external_gene_name), "loadings.cell.comp1"], row.names = colnames(pheat))
ann$Cluster <- ifelse(ann$Loadings > 0, "Recurrent", "Primary")
ann$Loadings <- abs(ann$Loadings) 
ann <- ann[order(stability.cell[match(colnames(pheat), stability.cell$external_gene_name), "loadings.cell.comp1"]),]
pheat <- pheat[, order(stability.cell[match(colnames(pheat), stability.cell$external_gene_name), "loadings.cell.comp1"])]
breaksList = seq(-10, 10, by = 0.5)
pheatmap(t(log2(pheat[6:10, ] / pheat[1:5, ])), scale = "none", cluster_rows =  F, cluster_cols = F, display_numbers = T, annotation_row = ann, , breaks = breaksList, legend_breaks = c(-8,-4,-1, 1, 4,8),  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)), gaps_row = table(ann$C)[1])
