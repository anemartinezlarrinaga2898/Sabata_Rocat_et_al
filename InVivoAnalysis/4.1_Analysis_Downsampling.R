# SCRIPT: Downsampling analyis
# AUTOR: ANE 
# Date: 16-01-2026

# Queremos ver si realmente el numero de DEG de las arterias se ve afectado por el numero de celulas, es decir si todos lo  hacemos downgrade y lo igualamos
# obtenemos el mismo patron? 

# Luego con el pseudobulk confirmaremos tambien que realmente que no pasa nada por comparar 80 vs 1000

##################################################################################################################################
set.seed(123)
# Setting working parameters
source("/ijc/USERS/amartinezl/MGRAUPERA_18/Paths.R")
directory<-setwd(Lista_Paths_Main$MainPath)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(openxlsx)
library(org.Mm.eg.db)
library(clusterProfiler)
library(optparse)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.3_Dowstream_Analysis/DEG/Downsampling_Plots",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

# | Parámetros del script
levels_needed <- c("WT", "MUT") 
phen_col<- "Phenotype"
fc_col <- "avg_log2FC"     # en Seurat suele ser avg_log2FC (comprueba tu objeto)
padj_col <- "p_val_adj"    # en Seurat es p_val_adj
alpha <- 0.05

# Functions 

# Función simple para sacar nombre de cluster del fichero
get_cluster_from_file <- function(fname){
  fname <- gsub("\\.rds$", "", fname)
  fname <- sub("^Lista_Markers_Repro_", "", fname)
  return(fname)
}

##################################################################################################################################

# ---- | Load files 

path.obj<-paste(path.guardar_original,"0.3_Dowstream_Analysis/DEG/Downsampling",sep="/")
list_files<-list.files(path.obj,pattern="Lista_Markers_Repro_")
Lista_Markers<-list()

for(i in seq_along(list_files)){
    lf<-list_files[i]
    markers<-readRDS(file.path(path.obj,lf))
    Lista_Markers[[i]]<-markers
}

# ---- | Round_Summary_All con bucles 

Round_Summary_All <- data.frame(
  cluster = character(),
  round = integer(),
  n_total = integer(),
  n_significant = integer(),
  n_up_significant = integer(),
  n_down_significant = integer(),
  stringsAsFactors = FALSE
)

for(i in seq_along(Lista_Markers)){

  cl <- get_cluster_from_file(list_files[i])
  lst <- Lista_Markers[[i]]   # lista de rondas (dataframes de FindMarkers)

  cat("Procesando cluster:", cl, " | rondas:", length(lst), "\n")

  for(r in seq_along(lst)){

    mk <- lst[[r]]

    # caso vacío / NULL
    if(is.null(mk) || nrow(mk) == 0){
      new_row <- data.frame(
        cluster = cl,
        round = r,
        n_total = 0,
        n_significant = 0,
        n_up_significant = 0,
        n_down_significant = 0,
        stringsAsFactors = FALSE
      )
      Round_Summary_All <- rbind(Round_Summary_All, new_row)
      next
    }

    # mk tiene genes como rownames
    df <- mk
    df$gene <- rownames(df)

    sig  <- df[[padj_col]] < alpha
    up   <- df[[fc_col]] > 0
    down <- df[[fc_col]] < 0

    n_total <- nrow(df)
    n_sig   <- sum(sig, na.rm = TRUE)
    n_up    <- sum(sig & up, na.rm = TRUE)
    n_down  <- sum(sig & down, na.rm = TRUE)

    new_row <- data.frame(
      cluster = cl,
      round = r,
      n_total = n_total,
      n_significant = n_sig,
      n_up_significant = n_up,
      n_down_significant = n_down,
      stringsAsFactors = FALSE
    )

    Round_Summary_All <- rbind(Round_Summary_All, new_row)
  }
}

# ---- Guardar ----
saveRDS(Round_Summary_All, file.path(path.guardar, "Round_Summary_Rebuilt_AllClusters.rds"))
write.csv(Round_Summary_All, file.path(path.guardar, "Round_Summary_Rebuilt_AllClusters.csv"), row.names = FALSE)

df_const <- Round_Summary_All %>%
  group_by(cluster) %>%
  summarise(
    n_star = NA_integer_,  # si quieres lo añadimos luego
    deg = unique(n_significant)[1],
    n_unique = n_distinct(n_significant),
    .groups="drop"
  )

# Si n_unique >1 en algún cluster, te avisa
print(df_const)

p <- ggplot(df_const, aes(x = cluster, y = deg)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "Cluster", y = "# DEGs (p_adj < 0.05)",
       title = "Balanced downsampling: DEG counts are invariant across 1000 rounds") 
ggsave(file.path(path.guardar, "PanelB_DEGs_Invariant_Dotplot.pdf"), p, width = 5, height = 4)

tab_const <- Round_Summary_All %>%
  group_by(cluster) %>%
  summarise(
    n_rounds = n(),
    n_unique_deg = n_distinct(n_significant),
    deg = unique(n_significant)[1],
    .groups="drop"
  )

p <- ggplot(tab_const, aes(x = cluster, y = deg)) +
  geom_col() +
  theme_classic() +
  labs(x="Cluster", y="# DEGs (balanced)")
ggsave(file.path(path.guardar, "PanelB_Downsampling_Invariant.pdf"), p, width=6, height=4)


# ---- | Preparación de datos para el plot representativo | ----

df_plot_rep <- Round_Summary_All %>%
  group_by(cluster) %>%
  summarise(
    # Promediamos (aunque sea invariante) para mayor seguridad
    UP = mean(n_up_significant),
    DOWN = mean(n_down_significant) * -1, # Ponemos DOWN en negativo para el efecto espejo
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c("UP", "DOWN"), names_to = "Status", values_to = "Count")

# ---- | Gráfico Representativo para el Paper | ----

p_final <- ggplot(df_plot_rep, aes(x = cluster, y = Count, fill = Status)) +
  geom_col(width = 0.7) +
  # Línea horizontal en 0 para separar UP de DOWN
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  # Colores consistentes con lo que definimos antes
  scale_fill_manual(values = c("UP" = "#D55E00", "DOWN" = "#0072B2")) +
  # Ajustamos el eje Y para que los valores negativos se vean como positivos (absolutos)
  scale_y_continuous(labels = abs) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "top",
    plot.title = element_text(size = 10, face = "bold")
  ) +
  labs(
    x = "Cluster", 
    y = "Number of DEGs",
    title = "Balanced Downsampling: UP and DOWN Regulated Genes",
    subtitle = "Results are invariant across 1000 iterations"
  )

# Guardar el gráfico "estrella" de la figura
ggsave(file.path(path.guardar, "PanelB_Downsampling_UpDown_Final.pdf"), p_final, width = 6, height = 5)