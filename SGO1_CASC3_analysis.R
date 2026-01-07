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
degs_raw <- read_excel("GSE153873/GSE153873_AD.vs.Old_diff.genes.xlsx")
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
# =====================================================
# ANOTACIÓN 450K - MÉTODO ESTABLE (SIN ERRORES S4)
# =====================================================

# 1. DMPs significativos
dmp_sig <- dmp_top[dmp_top$adj.P.Val < 0.05, ]
cat("✅ DMPs significativos:", nrow(dmp_sig), "\n")

# 2. ANOTACIÓN CORRECTA (método data.frame)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Cargar anotación como DATA.FRAME (evita error S4)
anno_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Convertir a data.frame legible
anno_df <- as.data.frame(anno_450k)
rownames(anno_df) <- anno_450k$Name

# Corregir nombres de probes 450K (puntos → guiones)
rownames(dmp_sig) <- gsub("\\.", "-", rownames(dmp_sig))

# Mapear SYMBOL
dmp_sig$SYMBOL <- anno_df[rownames(dmp_sig), "UCSC_RefGene_Name"]

# Limpiar SYMBOL
dmp_sig$SYMBOL <- sapply(strsplit(dmp_sig$SYMBOL, ";"), `[`, 1)
dmp_sig$SYMBOL[is.na(dmp_sig$SYMBOL) | dmp_sig$SYMBOL == ""] <- "NA"

cat("✅ SYMBOL asignados:", sum(dmp_sig$SYMBOL != "NA", na.rm=TRUE), "/", nrow(dmp_sig), "\n")

# 3. EXPORTAR Table2
table2_final <- dmp_sig[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "SYMBOL")]
colnames(table2_final)[1:6] <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
write.csv(table2_final, "Table2_DMPs_MTG.csv", row.names=TRUE)

# DEBUG ANOTACIÓN
head(rownames(dmp_sig))  # Probes ej "cg00000108"
head(rownames(anno_df))  # ¿Match?
sum(rownames(dmp_sig) %in% rownames(anno_df))  # Conteo match

# FIX 1: Fuerza match con Name col
dmp_sig$probeID <- rownames(dmp_sig)
anno_match <- anno_df[match(dmp_sig$probeID, anno_df$Name), "UCSC_RefGene_Name"]
dmp_sig$SYMBOL <- sapply(strsplit(anno_match, ";"), `[`, 1)
dmp_sig$SYMBOL[is.na(dmp_sig$SYMBOL) | dmp_sig$SYMBOL == "" | dmp_sig$SYMBOL == "c"] <- "NA"  # 'c' común en 450K

cat("FIXED SYMBOLs:", sum(dmp_sig$SYMBOL != "NA"), "/", nrow(dmp_sig), "\n")
cat("CASC3 probes:", sum(dmp_sig$SYMBOL == "CASC3"), "\n")

# Table2 con SYMBOL FIXED
table2_final <- dmp_sig[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "SYMBOL")]
write.csv(table2_final, "Table2_DMPs_MTG_fixed.csv", row.names=TRUE)
# =====================================================
# ANOTACIÓN 450K CORRECTA - PROBES YA ANOTADOS
# =====================================================

# 1. DMPs sig (ya con SYMBOL en rownames)
dmp_sig <- dmp_top[dmp_top$adj.P.Val < 0.05, ]
cat("✅ DMPs:", nrow(dmp_sig), "\n")

# 2. ANOTACIÓN: rownames(dmp_sig) = SYMBOL directo
dmp_sig$SYMBOL <- rownames(dmp_sig)  # ¡Ya son genes!
dmp_sig$probeID <- rownames(dmp_sig)  # Backup probe

# 3. Carga annotation para chr/pos (opcional, mejora Table3)
anno_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_df <- as.data.frame(anno_450k)

# Match SYMBOL a probe info (muchos probes → 1 gen)
probe_info <- anno_df[match(rownames(dmp_sig), anno_df$UCSC_RefGene_Name), 
                      c("chr", "pos", "UCSC_RefGene_Group")]

dmp_sig$chr <- probe_info$chr
dmp_sig$pos <- probe_info$pos  
dmp_sig$region <- probe_info$UCSC_RefGene_Group

cat("✅ Genes únicos:", length(unique(dmp_sig$SYMBOL[!is.na(dmp_sig$SYMBOL)])), "\n")

# 4. Table2 EXPORT (SYMBOL + chr/pos)
table2_final <- dmp_sig[, c("SYMBOL", "probeID", "chr", "pos", "region", 
                            "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
write.csv(table2_final, "Table2_DMPs_MTG.csv", row.names=FALSE)

head(table2_final)  # Ver HSPB3 etc.


# =========================================
# 5. GENE-LEVEL AGGREGATION (Table 3 FIJA)
# =========================================
dmp_genes <- dmp_sig %>%
  filter(!is.na(SYMBOL) & SYMBOL != "NA") %>%  
  group_by(SYMBOL) %>%
  summarise(
    nprobes = n(),
    meanlogFC = mean(logFC, na.rm=TRUE),
    minadjP = min(adj.P.Val, na.rm=TRUE),
    meanB = mean(B, na.rm=TRUE),
    .groups = 'drop'
  ) %>%
  arrange(minadjP) %>%
  filter(minadjP < 0.05)  # Solo FDR<0.05

write.csv(dmp_genes, "Table3_DMPs_gene_level.csv", row.names=FALSE)
casc3_row <- dmp_genes[dmp_genes$SYMBOL=="CASC3", ]
casc3_fc <- round(as.numeric(casc3_row$meanlogFC), 3)
casc3_fdr <- sprintf("%.3f", as.numeric(casc3_row$minadjP))
cat(sprintf("Manuscript: CASC3 MTG (logFC = %.2f, FDR = %.3f)", casc3_fc, as.numeric(casc3_fdr)))
cat("✅ Table3:", nrow(dmp_genes), "genes\n")
cat("🎉 CASC3: logFC =", casc3_fc, "| FDR =", casc3_fdr, "| nprobes =", casc3_row$nprobes, "\n")



# =========================================
# 6. SUPPLEMENTARY TABLES
# =========================================
# Supp Table S1: ALL probes FDR<0.05
write.csv(dmp_sig[, c("logFC","AveExpr","t","P.Value","adj.P.Val","B","SYMBOL")], 
          "Supp_Table_S1_DMPs_probes.csv", row.names=TRUE)

# Supp Table S2: ALL genes FDR<0.05
write.csv(dmp_genes, "Supp_Table_S2_DMPs_genes.csv", row.names=FALSE)

# =========================================
# 7. FIGURE 2 FINAL (con valores exactos)
# =========================================
casc3_mtg <- dmp_genes[dmp_genes$SYMBOL=="CASC3", ]
casc3_data <- data.frame(
  Tissue = c("Microglia", "MTG Cortex"),
  logFC = c(casc3_micro$`RA-RO`, casc3_mtg$meanlogFC),
  FDR = c(casc3_micro$FDR, casc3_mtg$minadjP)
)

p2 <- ggplot(casc3_data, aes(Tissue, logFC, fill=Tissue)) +
  geom_col(width=0.6, color="black") +
  geom_text(aes(label=sprintf("FDR=%.3f\nlogFC=%.2f", FDR, logFC)), 
            vjust=-0.5, size=4) +
  labs(title="CASC3 Replication: Microglia → MTG Cortex", 
       subtitle=paste("n_probes MTG =", casc3_mtg$nprobes),
       x="", y="logFC (AD-Control)") +
  theme_minimal() + theme(legend.position="none")

ggsave("Fig2_CASC3_replication.pdf", p2, width=8, height=5)

cat("🎉 COMPLETO! Push: git add .; git commit -m 'Table3 fixed + Supps'; git push\n")
