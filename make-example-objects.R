library(Seurat)
library(tidyverse)
library(zinbwave)
library(DESeq2)

load("data/seurat.rda")

## Make SE object from Seurat object
metadata = seurat@meta.data
metadata$genoregion = as.factor(paste(metadata$genotype, metadata$region, sep="_"))
counts = seurat@raw.data
counts = counts[Matrix::rowSums(counts >= 5) >= 5, ]
metadata$genoregion = gsub("anerior", "anterior", metadata$genoregion)
metadata$genotype = gsub("anerior", "anterior", metadata$genotype)
se = SummarizedExperiment(assays=list(counts=as.matrix(counts)),
                          colData=metadata)


## create and save dds object with ZINB-WaVE weights
weight_fn = file.path("data", "dds-weights.rda")
design = model.matrix(~genoregion, data=colData(se))
epsilon_de = 1e12
zinb = zinbFit(se, K=0, epsilon=epsilon_de, BPPARAM=BiocParallel::SerialParam(),
                X = design)
counts_zinb = zinbwave(se, fitted_model = zinb, K = 0, epsilon=epsilon_de)
dds_weights = DESeqDataSet(counts_zinb, design=~1+genoregion)
save(dds_weights, file=weight_fn)

## create and save dds object without weights
noweight_fn = file.path("data", "dds-noweights.rda")
ddsraw = DESeqDataSet(se, design=~1+genoregion)
save(ddsraw, file=noweight_fn)
