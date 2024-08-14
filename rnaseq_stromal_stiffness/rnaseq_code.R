suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(hgu95av2.db)
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(ExploreModelMatrix)
  library(cowplot)
  library(ComplexHeatmap)
  library(apeglm)
  library(gplots)
  library(ggplot2)
  library(microbenchmark)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(simplifyEnrichment)
})

counts <- read.table(file = "data/GSE179983_pdx_hs33_counts.txt", sep = ",", header = TRUE, row.names = 1)
coldata <- read.table("data/GSE179983_coldata.txt", sep = ",", header = TRUE, row.names = 1)

se <- SummarizedExperiment(assays = list(counts = as.matrix(counts)),colData = coldata)
se <- se[rowSums(assay(se, "counts")) > 5, ]
saveRDS(se, "data/GSE96870_se.rds")

dds <- DESeq2::DESeqDataSet(se,design = ~ stiffness)
dds <- DESeq(dds)

resStiffness <- results(dds, contrast = c("stiffness", "stiff", "soft"))
stiffnessDE <- as.data.frame(subset(resStiffness, padj < 0.05))
stiffnessDEgenes <- rownames(stiffnessDE)
stiffnessDEgenes <- sub('\\.[0-9]*$', '', stiffnessDEgenes)
gene_sets <- filter(msigdbr(species = "human"), gs_cat == "H" | gs_cat == "C6")

resStiffnessHallmark = enricher(gene = stiffnessDEgenes, 
                                TERM2GENE = gene_sets[, c("gs_name", "ensembl_gene")],
                                pvalueCutoff = 1,
                                qvalueCutoff = 1)
resStiffnessHallmarkTable = as.data.frame(resStiffnessHallmark)

barplot(resStiffnessHallmark, showCategory = 10)
ggplot(resStiffnessHallmarkTable[1:10, ],
       aes(x = -log(pvalue), y = factor(Description, levels = rev(Description)))) +
       geom_bar(stat = "identity") +
       geom_text(aes(x = pvalue, 
       label = sprintf("%.2e", p.adjust)), hjust = 1, col = "white") +
       ylab("")