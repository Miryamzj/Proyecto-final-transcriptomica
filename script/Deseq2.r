############################################################
# Análisis de expresión diferencial con DESeq2
# Proyecto: RNA-seq de hormigas infectadas por
# Ophiocordyceps unilateralis
############################################################

###############################
# Cargar librerías necesarias
###############################

library(DESeq2) # Expresión diferencial
library(tidyverse) # Manipulación de datos
library(pheatmap) # Heatmaps
library(apeglm) # Shrinkage de log2FC
library(EnhancedVolcano) # Volcano plots

###############################
# Importar matriz de conteos
###############################
# Leer tabla de conteos generada por featureCounts
counts <- read.delim("counts_Ant_clean.tsv", check.names = FALSE)
###############################
# Importar metadata
###############################

# Leer información experimental de las muestras
metadata <- read.csv("metadata_deseq2.csv", stringsAsFactors = FALSE)
###############################
# Preparar matriz de conteos
###############################

# Utilizar IDs génicos como nombres de fila
rownames(counts) <- counts$gene_id

# Eliminar columna gene_id
counts <- counts[, -1]

# Convertir a matriz numérica
counts <- as.matrix(counts)

# Asegurar que los conteos sean enteros
mode(counts) <- "integer"

###############################
# Preparar metadata
###############################

# Utilizar nombre de muestra como rowname
rownames(metadata) <- metadata$sample

# Definir condición experimental
metadata$condition <- factor(
  metadata$condition,
  levels = c("healthy", "infected")
)

# Reordenar metadata para coincidir con counts
metadata <- metadata[colnames(counts), ]

# Verificar coincidencia entre tablas
all(colnames(counts) == metadata$sample)

###############################
# Crear objeto DESeq2
###############################
# Diseño experimental:
# infected vs healthy

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~condition
)

###############################
# Filtrado de genes poco expresados
###############################
# Conservar genes con al menos
# 10 lecturas en 3 muestras

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("Genes retenidos:", nrow(dds), "\n")

###############################
# Ejecutar DESeq2
###############################
dds <- DESeq(dds)
###############################
# Obtener resultados
###############################

# Comparación:
# infected vs healthy

res <- results(
  dds,
  contrast = c(
    "condition",
    "infected",
    "healthy"
  )
)

###############################
# Shrinkage de log2FC
###############################
# Reduce ruido en genes poco expresados

res_shrink <- lfcShrink(
  dds,
  coef = "condition_infected_vs_healthy",
  type = "apeglm"
)

###############################
# Convertir a data.frame
###############################

res_df <- as.data.frame(res_shrink)
res_df$gene_id <- rownames(res_df)

###############################
# Ordenar por FDR
###############################

res_df <- res_df %>%
  arrange(padj)

###############################
# Definir umbrales
###############################
FDR <- 0.05
LFC <- 1

###############################
# Identificar DEGs
###############################

deg_sig <- res_df %>%
  filter(
    !is.na(padj),
    padj < FDR,
    abs(log2FoldChange) >= LFC
  )

cat(
  "Genes diferencialmente expresados:",
  nrow(deg_sig),
  "\n"
)

###############################
# Exportar tablas
###############################
write.csv(
  res_df,
  "DESeq2_all_genes.csv",
  row.names = FALSE
)

write.csv(
  deg_sig,
  "DESeq2_significant_genes.csv",
  row.names = FALSE
)

###############################
# Transformación VST
###############################
# Normalización para PCA y heatmaps

vsd <- vst(dds, blind = FALSE)

###############################
# PCA
###############################

vsd <- vst(dds, blind = FALSE) # Transformación VST para PCA

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE) # Extraer datos para PCA

percentVar <- round(100 * attr(pcaData, "percentVar")) # Porcentaje de varianza explicada por cada PC

library(ggplot2)
# Crear gráfico PCA con ggplot2
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition)) + # Colorear por condición
  geom_point(size = 4) + # Agregar puntos al gráfico
  scale_color_manual(
    # Definir colores para cada condición
    values = c(
      healthy = "blue", # Color para muestras sanas
      infected = "red" # Color para muestras infectadas
    )
  ) +
  xlab(
    paste0(
      "PC1: ", # Etiqueta para el eje X con porcentaje de varianza
      percentVar[1],
      "%"
    )
  ) +
  ylab(
    paste0(
      "PC2: ", # Etiqueta para el eje Y con porcentaje de varianza
      percentVar[2],
      "%"
    )
  ) +
  theme_bw()

jpeg(
  "PCA_condition.jpeg", # Guardar gráfico PCA como imagen JPEG
  width = 1800, # Ancho de la imagen en píxeles
  height = 1500, # Alto de la imagen en píxeles
  res = 300 # Resolución de la imagen en DPI
)

print(pca_plot) # Imprimir gráfico PCA en el dispositivo gráfico

dev.off() # Cerrar dispositivo gráfico para guardar la imagen

###############################
# MA Plot
###############################

jpeg(
  "MA_plot.jpg", # Guardar gráfico MA como imagen JPEG
  width = 2100, # Ancho de la imagen en píxeles
  height = 1800, # Alto de la imagen en píxeles
  res = 300 # Resolución de la imagen en DPI
)

# Definir colores
cols <- ifelse(
  # Genes significativos y log2FC > 0 en rojo, < 0 en azul, no significativos en gris
  res_shrink$padj < 0.05 & res_shrink$log2FoldChange > 0,
  "red",
  ifelse(
    res_shrink$padj < 0.05 & res_shrink$log2FoldChange < 0,
    "blue",
    "grey70"
  )
)

plotMA(
  # limites de los ejes, colores y título
  res_shrink,
  ylim = c(-8, 8),
  colSig = cols,
  cex = 0.7,
  main = "MA Plot"
)

legend(
  # Agregar leyendas para cada categoría de genes
  "topright",
  legend = c("UP", "DOWN", "Not significant"),
  col = c("red", "blue", "grey70"),
  pch = 16,
  bty = "n"
)

dev.off() # Cerrar dispositivo gráfico para guardar la imagen


###############################
# Volcano Plot
###############################

jpeg(
  # Guardar gráfico Volcano como imagen JPEG
  "Volcano_plot.jpg",
  width = 2400,
  height = 2000,
  res = 300
)

# Crear vector de colores por gen
keyvals <- ifelse(
  # Genes con log2FC <= -1 y FDR < 0.05 en azul, >= 1 y FDR < 0.05 en rojo, no significativos en gris
  res_shrink$log2FoldChange <= -1 & res_shrink$padj < 0.05,
  "blue",

  ifelse(
    # Genes con log2FC >= 1 y FDR < 0.05 en rojo, no significativos en gris
    res_shrink$log2FoldChange >= 1 & res_shrink$padj < 0.05,
    "red",
    "grey70"
  )
)

# Asignar nombres a la leyenda
names(keyvals)[keyvals == "red"] <- "UP"
names(keyvals)[keyvals == "blue"] <- "DOWN"
names(keyvals)[keyvals == "grey70"] <- "NO"

EnhancedVolcano(
  # Crear gráfico Volcano con EnhancedVolcano utilizando los resultados de DESeq2
  res_shrink,
  lab = NA,
  selectLab = NULL,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05, # Umbral de significancia para el eje Y (FDR)
  FCcutoff = 1, # Umbral de cambio de expresión para el eje X (log2FC)
  title = "Differential Expression Analysis",
  subtitle = "infected vs healthy",
  colCustom = keyvals,
  pointSize = 2.5, # Tamaño de los puntos en el gráfico
  drawConnectors = FALSE, # No dibujar líneas de conexión entre puntos y etiquetas
  boxedLabels = FALSE # No usar etiquetas con caja para los puntos seleccionados
)
dev.off()

###############################
# Heatmap de genes top
###############################
# librerias para una visualización más personalizada del heatmap
library(grid) # Para personalizar gráficos
library(RColorBrewer) # Para paletas de colores

metadata$label <- c(
  # Etiquetas para el heatmap
  "infected_43",
  "infected_44",
  "infected_45",
  "infected_46",
  "infected_47",
  "infected_48",
  "healthy_49",
  "healthy_50",
  "healthy_51",
  "healthy_52",
  "healthy_53",
  "healthy_54"
)

# Seleccionar top 50 genes
top50 <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 50)

# Extraer matriz transformada
mat <- assay(vsd)[
  top50$gene_id,
]

# 1. Asegurar nombres de columnas correctos
colnames(mat) <- metadata$label

# 2. Cerrar cualquier residuo gráfico en memoria
graphics.off()

# 3. Guardar heatmap con alta resolución
jpeg(
  "Heatmap_top50.jpg",
  width = 2400,
  height = 2200,
  res = 300
)

rownames(metadata) <- metadata$label # Asegurar que los nombres de fila de metadata coincidan con las columnas de mat
colnames(mat) <- metadata$label
attr(mat, "name") <- NULL

pheatmap(
  # Crear heatmap con pheatmap utilizando la matriz de expresión transformada para los top 50 genes
  mat, # Matriz de expresión transformada para los top 50 genes

  scale = "row", # Escalar por fila para visualizar patrones de expresión relativa

  fontsize_row = 8, # Tamaño de fuente para nombres de filas (genes)
  fontsize_col = 10, # Tamaño de fuente para nombres de columnas (muestras)

  clustering_distance_rows = "euclidean", # Distancia euclidiana para clustering de filas (genes)
  clustering_distance_cols = "euclidean", # Distancia euclidiana para clustering de columnas (muestras)
  clustering_method = "complete", # Método de clustering jerárquico completo

  color = colorRampPalette(
    # Paleta de colores personalizada para el heatmap
    rev(brewer.pal(11, "RdBu"))
  )(100),

  legend_breaks = c(-2, -1, 0, 1, 2), # Puntos de quiebre para la leyenda de colores
  name = "Z-score", # Etiqueta para la leyenda de colores
  legend_labels = c("-2", "-1", "0", "1", "2"), # Etiquetas para la leyenda de colores
  legend = TRUE, # Mostrar leyenda de colores

  main = "Top 50 Differentially Expressed Genes", # Título del heatmap

  angle_col = "90", # Rotar etiquetas de columnas 90 grados para mejor legibilidad

  border_color = NA # Sin bordes entre celdas para un aspecto más limpio
)

dev.off() # Cerrar dispositivo gráfico para guardar la imagen del heatmap

###############################
# Guardar información de sesión
###############################

writeLines(
  # Guardar información de sesión en un archivo de texto para reproducibilidad
  capture.output(
    sessionInfo()
  ),
  "sessionInfo.txt"
)
############################################################
# Fin del análisis
############################################################
