# =========================================
# SGO1/CASC3 Alzheimer's Disease Analysis Pipeline
# Reanalysis: GSE153873 (microglia RNA-seq) + GSE109887 (MTG methylome)
# Francisco Rodríguez-Esparragón, 2026 | medRxiv Supplementary Material
# =========================================

# 1. INSTALL AND LOAD REQUIRED PACKAGES
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("GEOquery","limma","minfi","IlluminaHumanMethylation450kanno.ilmn12.hg19"), update=FALSE)
install.packages(c("openxlsx","dplyr","ggplot2"), quiet=TRUE)

library(GEOquery)
library(openxlsx)
library(limma)
library(dplyr)
library(ggplot2)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# 2. GSE153873: MICROGLIA DIFFERENTIAL EXPRESSION (from precomputed Excel)
degs_raw <- read.xlsx("GSE153873_AD.vs.Old_diff.genes.xlsx", rowNames=FALSE)
colnames(degs_raw)[1] <- "symbol"
colnames(degs_raw)[6] <- "FDR"

degs_microglia <- degs_raw %>% filter(FDR < 0.05)
cat("✅ GSE153873: ", nrow(degs_microglia), "DEGs in AD microglia (FDR<0.05)\n")

# KEY RESULTS: SGO1 and CASC3
sgo1_micro <- degs_microglia[degs_microglia$symbol=="SGO1", c("symbol","RA-RO","PRA-RO","FDR")]
casc3_micro <- degs_microglia[degs_microglia$symbol=="CASC3", c("symbol","RA-RO","PRA-RO","FDR")]
print("SGO1 (microglia):"); print(sgo1_micro)
print("CASC3 (microglia):"); print(casc3_micro)

# EXPORT Table 1
write.csv(degs_microglia, "Table1_DEGs_microglia.csv", row.names=FALSE)

# 3. GSE109887: MTG DIFFERENTIAL METHYLATION ANALYSIS
gse109887 <- getGEO("GSE109887", GSEMatrix=TRUE, getGPL=FALSE)[[1]]
pheno <- pData(gse109887)
pheno$group <- ifelse(grepl("AD", pheno$title), "AD", "Control")

beta <- exprs(gse109887)
design <- model.matrix(~0 + pheno$group)
colnames(design) <- levels(factor(pheno$group))
cont_matrix <- makeContrasts(AD-Control, levels=design)

fit <- lmFit(beta, design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)
dmp_top <- topTable(fit2, coef=1, number=Inf, adjust.method="fdr")
dmp_sig <- dmp_top[dmp_top$adj.P.Val < 0.05, ]

cat("✅ GSE109887: ", nrow(dmp_sig), "significant DMPs in MTG cortex\n")
write.csv(dmp_sig, "Table2_DMPs_MTG.csv", row.names=TRUE)

# 4. GENE-LEVEL ANNOTATION (450K array)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probe_ann <- ann450k[rownames(dmp_sig), ]
dmp_sig$SYMBOL <- probe_ann$UCSC_RefGene_name

# Aggregate to gene level (mean logFC)
dmp_genes <- dmp_sig %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(
    mean_logFC = mean(logFC), 
    n_probes = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_logFC))

print("Top 10 MTG genes by |logFC|:")
print(head(dmp_genes[order(abs(dmp_genes$mean_logFC), decreasing=TRUE), ], 10))
print("CASC3 (MTG):"); print(dmp_genes[dmp_genes$SYMBOL=="CASC3", ])

write.csv(dmp_genes, "Table3_DMPs_gene_level.csv", row.names=FALSE)

# 5. FIGURE 1: VOLCANO PLOT (microglia DEGs)
volcano_data <- degs_raw %>%
  mutate(
    logFC = `RA-RO`, 
    neg_logP = -log10(`PRA-RO`),
    sig = case_when(
      FDR < 0.05 & abs(`RA-RO`) > 1 ~ "FDR<0.05 + |logFC|>1",
      FDR < 0.05 ~ "FDR<0.05", 
      TRUE ~ "NS"
    )
  )

p1 <- ggplot(volcano_data, aes(x=logFC, y=neg_logP, color=sig)) +
  geom_point(alpha=0.7, size=1) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue", alpha=0.7) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="blue", alpha=0.7) +
  scale_color_manual(values=c("FDR<0.05 + |logFC|>1"="red", 
                              "FDR<0.05"="orange", 
                              "NS"="grey")) +
  theme_minimal(base_size=12) +
  labs(
    title="GSE153873: AD Microglia DEGs",
    subtitle=paste0("SGO1 top-hit (logFC=", sgo1_micro$`RA-RO`, ") | n=", nrow(degs_microglia), " FDR<0.05"),
    x="log2 Fold Change (AD vs Aged Control)", 
    y="-log10(p-value)"
  ) +
  theme(legend.position="bottom", legend.title=element_blank())

ggsave("Fig1_Volcano_microglia.pdf", p1, width=9, height=6, dpi=300)
print(p1)

# 6. FIGURE 2: CASC3 REPLICATION BARPLOT
casc3_data <- data.frame(
  Tissue = c("Microglia (GSE153873)", "MTG Cortex (GSE109887)"),
  logFC = c(casc3_micro$`RA-RO`, dmp_genes[dmp_genes$SYMBOL=="CASC3", "mean_logFC"]),
  FDR = c(casc3_micro$FDR, dmp_genes[dmp_genes$SYMBOL=="CASC3", "mean_logFC"])
)

p2 <- ggplot(casc3_data, aes(x=Tissue, y=logFC, fill=Tissue)) +
  geom_col(alpha=0.8, color="black") +
  geom_text(aes(label=sprintf("FDR=%.3f", FDR)), vjust=-0.5) +
  theme_minimal(base_size=12) +
  labs(title="CASC3 Dysregulation: Microglia → Cortex Replication",
       x="", y="log2 Fold Change (AD vs Control)") +
  theme(legend.position="none")

ggsave("Fig2_CASC3_replication.pdf", p2, width=8, height=6, dpi=300)
print(p2)

# SUMMARY FOR METHODS SECTION
cat("\n🎉 ANALYSIS COMPLETE - medRxiv Supplementary Files Ready!\n")
cat("Generated files:\n")
print(list.files(pattern="^(Table|Fig).*"))
cat("\nMethods values for manuscript:\n")
cat("- Microglia DEGs (GSE153873):", nrow(degs_microglia), "\n")
cat("- SGO1:", sgo1_micro$symbol, "logFC =", sgo1_micro$`RA-RO`, "FDR =", sgo1_micro$FDR, "\n")
cat("- CASC3 (microglia): logFC =", casc3_micro$`RA-RO`, "FDR =", casc3_micro$FDR, "\n")
cat("- MTG DMPs (GSE109887):", nrow(dmp_sig), "\n")
cat("- CASC3 (MTG): logFC =", dmp_genes[dmp_genes$SYMBOL=="CASC3", "mean_logFC"], "\n")
