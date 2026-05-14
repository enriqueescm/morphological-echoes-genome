library(ggplot2)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(reshape2)

setwd("~/morphological-echoes-genome")

# Cargar datos
eigenvec <- read_table("data/processed/chr1_pca.eigenvec")
eigenval <- read_lines("data/processed/chr1_pca.eigenval") %>% as.numeric()
panel <- read_table("data/raw/integrated_call_samples_v3.20130502.ALL.panel")

# Varianza explicada
pve <- round(eigenval / sum(eigenval) * 100, 2)

# Unir con información de población
pca_data <- eigenvec %>%
  rename(sample = `#IID`) %>%
  left_join(panel, by = "sample")

# Colores por superpoblación
colores <- c(
  "AFR" = "#E41A1C",
  "EUR" = "#377EB8",
  "EAS" = "#4DAF4A",
  "SAS" = "#FF7F00",
  "AMR" = "#984EA3"
)

# Colores para 26 poblaciones
colores_pop <- c(
  "YRI" = "#E41A1C", "LWK" = "#FF4444", "GWD" = "#CC0000",
  "MSL" = "#FF6666", "ESN" = "#AA0000", "ASW" = "#FF9999",
  "ACB" = "#FFCCCC",
  "CEU" = "#377EB8", "TSI" = "#4444FF", "FIN" = "#0000CC",
  "GBR" = "#6666FF", "IBS" = "#99AAFF",
  "CHB" = "#4DAF4A", "JPT" = "#00CC00", "CHS" = "#009900",
  "CDX" = "#66FF66", "KHV" = "#336633",
  "GIH" = "#FF7F00", "PJL" = "#FF9933", "BEB" = "#CC6600",
  "STU" = "#FFAA44", "ITU" = "#FFD700",
  "MXL" = "#984EA3", "PUR" = "#CC66FF", "CLM" = "#7700AA",
  "PEL" = "#DD99FF"
)

# Centroides por superpoblación
centroides <- pca_data %>%
  group_by(super_pop) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

# Plot 1 — superpoblaciones
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = super_pop, fill = super_pop)) +
  geom_point(alpha = 0.7, size = 2, shape = 21, stroke = 0.3) +
  scale_color_manual(values = colores, name = "Superpopulation") +
  scale_fill_manual(values = colores, name = "Superpopulation") +
  labs(
    title = "Population Structure — 1000 Genomes Project",
    subtitle = paste0("Chr1 | ", nrow(pca_data), " individuals | 22,415 SNPs after LD pruning"),
    x = paste0("PC1 (", pve[1], "% variance explained)"),
    y = paste0("PC2 (", pve[2], "% variance explained)"),
    caption = "Source: 1000 Genomes Project Phase 3"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p1
ggsave("figures/01_PCA_superpopulations.png", p1, width = 10, height = 7, dpi = 300)

# Plot 2 — 26 poblaciones
p2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = pop, fill = pop)) +
  geom_point(alpha = 0.7, size = 2, shape = 21, stroke = 0.3) +
  scale_color_manual(values = colores_pop, name = "Population") +
  scale_fill_manual(values = colores_pop, name = "Population") +
  labs(
    title = "Population Structure — 26 Populations",
    subtitle = paste0("Chr1 | ", nrow(pca_data), " individuals | 22,415 SNPs after LD pruning"),
    x = paste0("PC1 (", pve[1], "% variance explained)"),
    y = paste0("PC2 (", pve[2], "% variance explained)"),
    caption = "Source: 1000 Genomes Project Phase 3"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p2
ggsave("figures/02_PCA_populations.png", p2, width = 12, height = 7, dpi = 300)

# Plot 3 — etiquetas por superpoblación
p3 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = super_pop, fill = super_pop)) +
  geom_point(alpha = 0.5, size = 1.8, shape = 21, stroke = 0.2) +
  geom_text(data = centroides, aes(label = super_pop),
            size = 4.5, fontface = "bold",
            color = "black",
            show.legend = FALSE) +
  scale_color_manual(values = colores, name = "Superpopulation") +
  scale_fill_manual(values = colores, name = "Superpopulation") +
  labs(
    title = "Human Genetic Variation Mirrors Geography",
    subtitle = paste0("Chr1 | ", nrow(pca_data), " individuals | 22,415 SNPs after LD pruning"),
    x = paste0("PC1 (", pve[1], "% variance explained)"),
    y = paste0("PC2 (", pve[2], "% variance explained)"),
    caption = "Source: 1000 Genomes Project Phase 3"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p3
ggsave("figures/03_PCA_labeled.png", p3, width = 10, height = 7, dpi = 300)

# Plot 4 — FST heatmap
fst_data <- read_table("data/processed/chr1_fst.fst.summary")

pops <- unique(c(fst_data$`#POP1`, fst_data$POP2))
fst_matrix <- matrix(0, nrow = length(pops), ncol = length(pops))
rownames(fst_matrix) <- pops
colnames(fst_matrix) <- pops

for (i in 1:nrow(fst_data)) {
  p1 <- fst_data$`#POP1`[i]
  p2 <- fst_data$POP2[i]
  val <- fst_data$HUDSON_FST[i]
  fst_matrix[p1, p2] <- val
  fst_matrix[p2, p1] <- val
}

fst_long <- melt(fst_matrix)
colnames(fst_long) <- c("Pop1", "Pop2", "FST")

p4 <- ggplot(fst_long, aes(x = Pop1, y = Pop2, fill = FST)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(FST, 3)), size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "#FFFFD4", high = "#E41A1C",
                      name = "FST",
                      limits = c(0, 0.16)) +
  labs(
    title = "Genetic Differentiation Between Human Populations",
    subtitle = "Hudson FST — Chr1 | 22,415 SNPs | 1000 Genomes Project Phase 3",
    x = NULL, y = NULL,
    caption = "Higher FST = greater genetic differentiation"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    panel.grid = element_blank()
  )

p4
ggsave("figures/04_FST_heatmap.png", p4, width = 8, height = 7, dpi = 300)

# Análisis de genes candidatos craneofaciales en Chr1
# RUNX2: chr1:45,392,147-45,746,609
# COL1A1: chr1:101,941,859-101,974,317

# Cargar FST por SNP AFR-EAS (par más diferenciado)
fst_snp <- read_table("data/processed/chr1_fst_snp.AFR.EAS.fst.var")

# Definir regiones de genes candidatos
genes_candidatos <- data.frame(
  gene = c("RUNX2", "COL1A1"),
  start = c(45392147, 101941859),
  end = c(45746609, 101974317)
)

# Etiquetar SNPs en genes candidatos
fst_snp <- fst_snp %>%
  mutate(gene = case_when(
    POS >= 45392147 & POS <= 45746609 ~ "RUNX2",
    POS >= 101941859 & POS <= 101974317 ~ "COL1A1",
    TRUE ~ NA_character_
  ))

# Manhattan-style plot
p5 <- ggplot(fst_snp, aes(x = POS/1e6, y = HUDSON_FST)) +
  geom_point(alpha = 0.3, size = 0.8, color = "grey60") +
  geom_point(data = filter(fst_snp, !is.na(gene)),
             aes(color = gene), size = 2.5, alpha = 0.9) +
  geom_hline(yintercept = quantile(fst_snp$HUDSON_FST, 0.99, na.rm = TRUE),
             linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_color_manual(values = c("RUNX2" = "#E41A1C", "COL1A1" = "#377EB8"),
                     name = "Candidate Gene") +
  labs(
    title = "FST Along Chromosome 1 — AFR vs EAS",
    subtitle = "Craniofacial candidate genes highlighted | 1000 Genomes Project Phase 3",
    x = "Position (Mb)",
    y = "Hudson FST",
    caption = "Dashed line = 99th percentile FST threshold"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p5
ggsave("figures/05_FST_candidates_chr1.png", p5, width = 12, height = 6, dpi = 300)