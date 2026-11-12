# =========================================================
# STATIS DUAL NEUTROSOFICO PARA COOPERATIVAS DE AHORRO Y CREDITO
# Script adaptado de la metodologia STATIS Dual Neutrosófico
# Datos esperados: DATOS_GABI_TUCKER3.xlsx
# =========================================================
# Este script implementa la metodología STATIS Dual Neutrosófica:
# Fase 1: Preparación de datos neutrosóficos (mismo que STATIS clásico)
# Fase 2: Inter-estructura Dual (matrices entre variables, no entre individuos)
# Fase 3: Compromiso Dual Neutrosófico
# Fase 4: Intra-estructura Dual (ACP del compromiso entre variables)
# Fase 5: Trayectorias de la estructura de variables en el tiempo
# =========================================================

# ------------------------------ 0) PAQUETES ------------------------------
pkgs <- c(
  "readxl", "dplyr", "tidyr", "stringr", "ggplot2", "ggrepel",
  "purrr", "tibble", "reshape2"
)
new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs) > 0) install.packages(new_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
set.seed(123)

# ------------------------------ 1) PARAMETROS GENERALES ------------------------------
ruta_archivo <- "DATOS_GABI_TUCKER3.xlsx"
vars <- c("DEP", "SCC", "INTFIN", "LIQCOR", "MORAMP", "ROA", "ROE")
years <- 2016:2023

# Pesos para deneutrosificacion (mismos que en STATIS clásico)
wT <- 0.50
wI <- 0.30
wF <- 0.20

# Numero de dimensiones principales del espacio de compromiso dual
nf <- 3

# ------------------------------ 2) FUNCIONES AUXILIARES ------------------------------
rescale01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (any(is.infinite(rng)) || is.na(diff(rng)) || diff(rng) == 0) {
    return(rep(0.5, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

safe_sd <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) 1 else s
}

calc_target_truth <- function(x, bench) {
  dev <- abs(x - bench)
  max_dev <- max(dev, na.rm = TRUE)
  if (is.na(max_dev) || max_dev == 0) return(rep(1, length(x)))
  1 - dev / max_dev
}

calc_optimal_truth <- function(x, target = 1) {
  dev <- abs(x - target)
  max_dev <- max(dev, na.rm = TRUE)
  if (is.na(max_dev) || max_dev == 0) return(rep(1, length(x)))
  1 - dev / max_dev
}

center_scale_matrix <- function(X) {
  X <- as.matrix(X)
  Xs <- scale(X, center = TRUE, scale = TRUE)
  Xs[is.na(Xs)] <- 0
  Xs
}

# Coeficiente RV (mismo que en STATIS clásico)
rv_coeff <- function(S1, S2) {
  num <- sum(S1 * S2)
  den <- sqrt(sum(S1 * S1) * sum(S2 * S2))
  if (is.na(den) || den == 0) return(0)
  as.numeric(num / den)
}

# Funcion para construir matriz anual (cooperativas x variables)
make_year_matrix <- function(df_long, year_value, value_col, coops, vars_ord) {
  tmp <- df_long %>%
    filter(year == year_value) %>%
    select(coop, variable, !!sym(value_col)) %>%
    tidyr::pivot_wider(names_from = variable, values_from = !!sym(value_col)) %>%
    arrange(match(coop, coops))
  
  X <- as.matrix(tmp[, vars_ord, drop = FALSE])
  rownames(X) <- tmp$coop
  colnames(X) <- vars_ord
  X
}

# ------------------------------ 3) LECTURA Y REESTRUCTURACION ------------------------------
col_names <- c("coop", unlist(lapply(years, function(y) paste(vars, y, sep = "_"))))
raw <- read_excel(ruta_archivo, col_names = FALSE)
names(raw) <- col_names
raw <- raw[-1, ]
raw$coop <- as.character(raw$coop)
raw[-1] <- lapply(raw[-1], as.numeric)

long <- raw %>%
  pivot_longer(-coop, names_to = "var_year", values_to = "value") %>%
  tidyr::extract(var_year, into = c("variable", "year"), regex = "(.*)_(\\d+)", convert = TRUE) %>%
  arrange(coop, variable, year)

coops <- unique(long$coop)

cat("====================================\n")
cat("ESTRUCTURA DE DATOS PARA STATIS DUAL NEUTROSOFICO\n")
cat("====================================\n")
cat("Numero de cooperativas:", length(coops), "\n")
cat("Numero de variables:", length(vars), "\n")
cat("Numero de anios:", length(years), "\n\n")

# ------------------------------ 4) CONSTRUCCION DE T, I, F (igual que STATIS clásico) ------------------------------
long <- long %>%
  group_by(variable, year) %>%
  mutate(
    norm = rescale01(value),
    bench_med = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(coop, variable) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    d_prev = abs(value - lag(value)),
    d_next = abs(lead(value) - value),
    vol_local_raw = rowMeans(cbind(ifelse(is.na(d_prev), 0, d_prev),
                                   ifelse(is.na(d_next), 0, d_next))),
    vol_local = rescale01(vol_local_raw)
  ) %>%
  ungroup() %>%
  group_by(variable, year) %>%
  mutate(
    dev_sector_raw = abs(value - bench_med) / safe_sd(value),
    dev_sector = rescale01(dev_sector_raw)
  ) %>%
  ungroup()

long <- long %>%
  mutate(
    T = case_when(
      variable %in% c("ROA", "ROE") ~ norm,
      variable == "MORAMP" ~ 1 - norm,
      variable == "LIQCOR" ~ calc_target_truth(value, bench_med),
      variable %in% c("DEP", "SCC") ~ pmin(1, pmax(0, (value - q25) / pmax(q75 - q25, 1e-8))),
      variable == "INTFIN" ~ calc_optimal_truth(value, target = 1.0),
      TRUE ~ norm
    ),
    I = case_when(
      variable %in% c("ROA", "ROE") ~ pmin(1, 0.70 * vol_local + 0.30 * dev_sector),
      variable == "MORAMP" ~ pmin(1, 0.60 * vol_local + 0.40 * dev_sector),
      variable == "LIQCOR" ~ pmin(1, 0.50 * vol_local + 0.50 * dev_sector),
      variable %in% c("DEP", "SCC") ~ pmin(1, 0.65 * vol_local + 0.35 * dev_sector),
      variable == "INTFIN" ~ pmin(1, 0.60 * vol_local + 0.40 * dev_sector),
      TRUE ~ pmin(1, 0.70 * vol_local + 0.30 * dev_sector)
    )
  ) %>%
  mutate(
    T = pmin(pmax(T, 0), 1),
    I = pmin(pmax(I, 0), 1),
    F = pmin(pmax(1 - T, 0), 1),
    score_neutro = wT * T + wI * (1 - I) + wF * (1 - F)
  )

# ------------------------------ 5) MATRICES ANUALES PARA STATIS DUAL ------------------------------
X_years <- lapply(years, function(yy) make_year_matrix(long, yy, "score_neutro", coops, vars))
names(X_years) <- as.character(years)

# Centrar y escalar cada tabla
X_years_cs <- lapply(X_years, center_scale_matrix)

# ------------------------------ 6) INTER-ESTRUCTURA DUAL: MATRICES ENTRE VARIABLES ------------------------------
# DIFERENCIA CLAVE CON STATIS CLASICO: W_k = t(X_k) %*% X_k (variables x variables)
W_k_dual <- lapply(X_years_cs, function(X) t(X) %*% X)
names(W_k_dual) <- as.character(years)

# ------------------------------ 7) COEFICIENTE RV DUAL ------------------------------
# Matriz de similitud entre años basada en la estructura de variables
K <- length(years)
RV_dual <- matrix(0, nrow = K, ncol = K, dimnames = list(as.character(years), as.character(years)))

for (i in 1:K) {
  for (j in 1:K) {
    RV_dual[i, j] <- rv_coeff(W_k_dual[[i]], W_k_dual[[j]])
  }
}

cat("\n====================================\n")
cat("MATRIZ RV DUAL (Similitud entre estructuras de variables)\n")
cat("====================================\n")
print(round(RV_dual, 6))

# ------------------------------ 8) PESOS DUALES (COMPROMISO) ------------------------------
eig_RV_dual <- eigen(RV_dual)
weights_dual <- eig_RV_dual$vectors[, 1]
weights_dual <- as.numeric(weights_dual / sum(weights_dual))
if (sum(weights_dual < 0) > length(weights_dual) / 2) {
  weights_dual <- -weights_dual
  weights_dual <- weights_dual / sum(weights_dual)
}
names(weights_dual) <- as.character(years)

cat("\n====================================\n")
cat("PESOS STATIS DUAL POR ANIO\n")
cat("====================================\n")
print(round(data.frame(year = years, weight = weights_dual), 6))

# ------------------------------ 9) MATRIZ DE COMPROMISO DUAL ------------------------------
# Compromiso: combinación convexa de las matrices de variables
W_compromise_dual <- Reduce(`+`, Map(function(Wk, a) a * Wk, W_k_dual, weights_dual))
rownames(W_compromise_dual) <- vars
colnames(W_compromise_dual) <- vars

cat("\n====================================\n")
cat("MATRIZ DE COMPROMISO DUAL (Estructura consensuada de variables)\n")
cat("====================================\n")
print(round(W_compromise_dual, 6))

# ------------------------------ 10) INTRA-ESTRUCTURA DUAL: ACP DEL COMPROMISO ------------------------------
eig_comp_dual <- eigen(W_compromise_dual)

# Seleccionar dimensiones con autovalor positivo
positive_idx_dual <- which(eig_comp_dual$values > 1e-8)
if (length(positive_idx_dual) < nf) {
  nf <- length(positive_idx_dual)
}

lambda_dual <- eig_comp_dual$values[positive_idx_dual][1:nf]
U_dual <- eig_comp_dual$vectors[, positive_idx_dual, drop = FALSE][, 1:nf, drop = FALSE]
rownames(U_dual) <- vars
colnames(U_dual) <- paste0("Dim", 1:nf)

# Coordenadas de las variables en el espacio de compromiso dual
coord_vars_dual <- as.data.frame(U_dual %*% diag(sqrt(lambda_dual), nrow = length(lambda_dual)))
colnames(coord_vars_dual) <- paste0("Dim", 1:nf)
coord_vars_dual$variable <- vars
coord_vars_dual$inercia_contrib <- rowSums(coord_vars_dual[, paste0("Dim", 1:nf), drop = FALSE]^2)

# Porcentajes de inercia explicada
inercia_pct_dual <- 100 * lambda_dual / sum(eig_comp_dual$values[eig_comp_dual$values > 0])

cat("\n====================================\n")
cat("INERCIA DEL ESPACIO DE COMPROMISO DUAL\n")
cat("====================================\n")
print(data.frame(
  dimension = paste0("Dim", 1:nf),
  eigenvalue = round(lambda_dual, 6),
  inertia_pct = round(inercia_pct_dual, 6)
))

# ------------------------------ 11) PROYECCION DE AÑOS INDIVIDUALES (TRAYECTORIAS DE LA ESTRUCTURA) ------------------------------
# Proyectar cada matriz W_k en el espacio de compromiso
partial_dual_coords <- lapply(seq_along(years), function(i) {
  # F_k = W_k * U * Lambda^{-1/2}
  inv_sqrt_lambda <- diag(1 / sqrt(lambda_dual), nrow = length(lambda_dual))
  Fk <- W_k_dual[[i]] %*% U_dual %*% inv_sqrt_lambda
  out <- as.data.frame(Fk)
  colnames(out) <- paste0("Dim", 1:nf)
  out$variable <- vars
  out$year <- years[i]
  out$weight <- weights_dual[i]
  out
})

partial_dual_df <- bind_rows(partial_dual_coords)

# ------------------------------ 12) COORDENADAS DE COOPERATIVAS COMO SUPLEMENTARIAS ------------------------------
# Proyectar las cooperativas en el espacio de compromiso dual (elementos suplementarios)
# Para cada año, calcular coordenadas de cooperativas: G_k = X_k * U * Lambda^{-1/2}
coops_sup_list <- lapply(seq_along(years), function(i) {
  inv_sqrt_lambda <- diag(1 / sqrt(lambda_dual), nrow = length(lambda_dual))
  Gk <- X_years_cs[[i]] %*% U_dual %*% inv_sqrt_lambda
  out <- as.data.frame(Gk)
  colnames(out) <- paste0("Dim", 1:nf)
  out$coop <- coops
  out$year <- years[i]
  out
})

coops_sup_df <- bind_rows(coops_sup_list)

# Coordenadas promedio de cooperativas en el compromiso dual (consenso)
coops_consensus <- coops_sup_df %>%
  group_by(coop) %>%
  summarise(
    Dim1 = mean(Dim1, na.rm = TRUE),
    Dim2 = mean(Dim2, na.rm = TRUE),
    Dim3 = mean(Dim3, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------ 13) CLUSTER DE VARIABLES ------------------------------
set.seed(123)
k_clusters_vars <- 3  # Número de clusters para variables
km_vars <- kmeans(coord_vars_dual[, paste0("Dim", 1:min(2, nf))], centers = k_clusters_vars, nstart = 25)
coord_vars_dual$cluster <- paste0("Cluster_", km_vars$cluster)

# -----------------------------------------------------------------------------
# AJUSTE EXCLUSIVAMENTE GRAFICO PARA LOS CIRCULOS DE CORRELACION
# -----------------------------------------------------------------------------
# Las matrices W_k se calculan como Z_k^T Z_k. Debido a que scale() utiliza la
# desviacion estandar muestral, con n cooperativas la diagonal de W_k es n - 1.
# Por ello, las coordenadas U Lambda^(1/2) se reescalan SOLO PARA DIBUJAR los
# mapas con circulo unitario. Ningun calculo STATIS, RV, peso, autovalor,
# contribucion, calidad de representacion o trayectoria es modificado.

graph_scale_factor <- sqrt(length(coops) - 1)

coord_vars_dual_plot <- coord_vars_dual %>%
  mutate(
    Dim1_plot = Dim1 / graph_scale_factor,
    Dim2_plot = Dim2 / graph_scale_factor,
    Dim3_plot = Dim3 / graph_scale_factor
  )

cat("Factor de reescalado grafico para circulo unitario:",
    graph_scale_factor, "\n")

# ------------------------------ 14) EXPORTAR TABLAS ------------------------------
write.csv(long, "datos_neutrosoficos_statis_dual_long.csv", row.names = FALSE)
write.csv(data.frame(year = years, weight = weights_dual), "pesos_statis_dual_neutrosofico.csv", row.names = FALSE)
write.csv(as.data.frame(RV_dual) %>% mutate(year = rownames(RV_dual)) %>% relocate(year), "matriz_rv_dual_neutrosofica.csv", row.names = FALSE)
write.csv(coord_vars_dual, "coordenadas_variables_compromiso_dual.csv", row.names = FALSE)
write.csv(partial_dual_df, "trayectorias_estructura_variables_dual.csv", row.names = FALSE)
write.csv(coops_consensus, "coordenadas_cooperativas_compromiso_dual.csv", row.names = FALSE)
write.csv(data.frame(dimension = paste0("Dim", 1:nf), eigenvalue = lambda_dual, inertia_pct = inercia_pct_dual),
          "inercia_compromiso_dual.csv", row.names = FALSE)

# ------------------------------ 15) GRAFICOS ------------------------------

# 15.1 Heatmap de similitud RV Dual
RV_dual_df <- as.data.frame(as.table(RV_dual))
colnames(RV_dual_df) <- c("year1", "year2", "RV")

p_rv_dual <- ggplot(RV_dual_df, aes(x = year1, y = year2, fill = RV)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", RV)), size = 3) +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.4) +
  labs(title = "Matriz de similitud RV Dual Neutrosófico",
       subtitle = "Similitud entre estructuras de relaciones de variables financieras",
       x = "Año", y = "Año", fill = "RV") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("heatmap_rv_dual_neutrosofico.png", p_rv_dual, width = 8, height = 7, dpi = 300)

# 15.2 Mapa de compromiso dual (circulo de correlaciones de variables)
# CORRECCION: las coordenadas se dividen por sqrt(n - 1) SOLO para dibujarlas
# en la misma escala que el circulo unitario.
p_vars_dual <- ggplot(
  coord_vars_dual_plot,
  aes(x = Dim1_plot, y = Dim2_plot, color = cluster)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  annotate(
    "path",
    x = cos(seq(0, 2 * pi, length.out = 200)),
    y = sin(seq(0, 2 * pi, length.out = 200)),
    color = "grey50",
    linetype = "dotted"
  ) +
  geom_segment(
    aes(x = 0, y = 0, xend = Dim1_plot, yend = Dim2_plot),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 0.8,
    color = "darkblue"
  ) +
  geom_point(size = 4) +
  # Etiquetas generales: se mantienen como en la figura original
  ggrepel::geom_text_repel(
    data = coord_vars_dual_plot %>%
      filter(!variable %in% c("INTFIN", "LIQCOR")),
    aes(label = variable),
    size = 5,
    fontface = "bold",
    max.overlaps = 100
  ) +
  # Correccion solicitada por la tutora:
  # INTFIN y LIQCOR se colocan centradas y por debajo de sus puntos
  geom_text(
    data = coord_vars_dual_plot %>%
      filter(variable %in% c("INTFIN", "LIQCOR")),
    aes(label = variable),
    nudge_y = -0.07,
    vjust = 1,
    size = 5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  labs(
    title = "Neutrosophic-Score Dual STATIS Compromise Map",
    subtitle = sprintf(
      paste0(
        "Consensus structure of relationships among financial variables\n",
        "Inertia: Dim1 = %.2f%%, Dim2 = %.2f%%"
      ),
      inercia_pct_dual[1],
      inercia_pct_dual[2]
    ),
    x = paste0("Dimension 1 (", round(inercia_pct_dual[1], 2), "%)"),
    y = paste0("Dimension 2 (", round(inercia_pct_dual[2], 2), "%)"),
    color = "Cluster"
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(-1.05, 1.05),
    ylim = c(-1.05, 1.05)
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_vars_dual)

ggsave(
  "mapa_compromiso_dual_variables.png",
  p_vars_dual,
  width = 10, height = 8, dpi = 300
)

ggsave(
  "Figure_1_Neutrosophic_Dim1_Dim2_CORRECTED.pdf",
  p_vars_dual,
  device = cairo_pdf,
  width = 10, height = 8, units = "in"
)

# 15.3 Trayectorias de la estructura de variables (evolución temporal)
p_traj_estructura <- ggplot(partial_dual_df, aes(x = Dim1, y = Dim2, group = variable, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_path(linewidth = 0.8, alpha = 0.7) +
  geom_point(aes(shape = factor(year)), size = 2.5) +
  ggrepel::geom_text_repel(data = partial_dual_df %>% filter(year == max(year)),
                           aes(label = variable), size = 4, fontface = "bold", max.overlaps = 100) +
  labs(title = "Trayectorias de la Estructura de Variables Financieras (2016-2023)",
       subtitle = "Evolución de la posición de cada variable en el espacio de compromiso dual",
       x = paste0("Dimensión 1 (", round(inercia_pct_dual[1], 2), "%)"),
       y = paste0("Dimensión 2 (", round(inercia_pct_dual[2], 2), "%)"),
       shape = "Año") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("trayectorias_estructura_variables_dual.png", p_traj_estructura, width = 11, height = 8, dpi = 300)

# 15.4 Biplot dual: variables + cooperativas (suplementarias)
# Tomar solo las cooperativas más representativas para no saturar el gráfico
coops_contrib <- coops_consensus %>%
  mutate(contrib = sqrt(Dim1^2 + Dim2^2)) %>%
  arrange(desc(contrib)) %>%
  slice(1:12)

p_biplot_dual <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  # Variables activas (vectores)
  geom_segment(data = coord_vars_dual, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.9, color = "darkred") +
  geom_point(data = coord_vars_dual, aes(x = Dim1, y = Dim2), size = 4, color = "darkred", shape = 17) +
  ggrepel::geom_text_repel(data = coord_vars_dual, aes(x = Dim1, y = Dim2, label = variable),
                           size = 5, fontface = "bold", color = "darkred", max.overlaps = 100) +
  # Cooperativas suplementarias
  geom_point(data = coops_contrib, aes(x = Dim1, y = Dim2), size = 2.5, alpha = 0.7, color = "steelblue") +
  ggrepel::geom_text_repel(data = coops_contrib, aes(x = Dim1, y = Dim2, label = coop),
                           size = 3, color = "steelblue", max.overlaps = 100) +
  labs(title = "Biplot Dual Neutrosófico",
       subtitle = "Variables activas (vectores rojos) y cooperativas suplementarias (puntos azules)",
       x = paste0("Dimensión 1 (", round(inercia_pct_dual[1], 2), "%)"),
       y = paste0("Dimensión 2 (", round(inercia_pct_dual[2], 2), "%)")) +
  coord_fixed(ratio = 1) +
  theme_minimal()

ggsave("biplot_dual_variables_cooperativas.png", p_biplot_dual, width = 11, height = 8, dpi = 300)

# 15.5 Contribución de variables por dimensión
contrib_vars <- coord_vars_dual %>%
  mutate(
    contrib_dim1 = abs(Dim1) / sum(abs(Dim1)) * 100,
    contrib_dim2 = abs(Dim2) / sum(abs(Dim2)) * 100
  ) %>%
  select(variable, contrib_dim1, contrib_dim2, cluster) %>%
  arrange(desc(contrib_dim1))

p_contrib <- ggplot(contrib_vars, aes(x = reorder(variable, contrib_dim1), y = contrib_dim1, fill = cluster)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Contribución de Variables a la Dimensión 1 del Compromiso Dual",
       x = "Variable", y = "Contribución (%)") +
  theme_minimal()

ggsave("contribucion_variables_dim1_dual.png", p_contrib, width = 8, height = 6, dpi = 300)

# 15.6 Evolución de la correlación entre variables clave en el tiempo
# Extraer elementos de W_k_dual (correlaciones entre variables por año)
cor_evolution <- lapply(seq_along(years), function(i) {
  W <- W_k_dual[[i]]
  # Normalizar a matriz de correlación
  d <- sqrt(diag(W))
  cor_matrix <- W / outer(d, d)
  data.frame(
    year = years[i],
    var1 = rep(rownames(cor_matrix), each = ncol(cor_matrix)),
    var2 = rep(colnames(cor_matrix), times = nrow(cor_matrix)),
    cor = as.vector(cor_matrix)
  )
}) %>% bind_rows()

# Filtrar pares relevantes (ej: ROA-ROE, ROA-MORAMP, DEP-SCC, INTFIN-ROA)
pairs_of_interest <- c("ROA_ROE", "ROA_MORAMP", "DEP_SCC", "INTFIN_ROA", "INTFIN_ROE")
cor_evolution_filtered <- cor_evolution %>%
  filter(paste(var1, var2, sep = "_") %in% pairs_of_interest | 
           paste(var2, var1, sep = "_") %in% pairs_of_interest) %>%
  mutate(pair = case_when(
    (var1 == "ROA" & var2 == "ROE") | (var1 == "ROE" & var2 == "ROA") ~ "ROA-ROE",
    (var1 == "ROA" & var2 == "MORAMP") | (var1 == "MORAMP" & var2 == "ROA") ~ "ROA-MORAMP",
    (var1 == "DEP" & var2 == "SCC") | (var1 == "SCC" & var2 == "DEP") ~ "DEP-SCC",
    (var1 == "INTFIN" & var2 == "ROA") | (var1 == "ROA" & var2 == "INTFIN") ~ "INTFIN-ROA",
    (var1 == "INTFIN" & var2 == "ROE") | (var1 == "ROE" & var2 == "INTFIN") ~ "INTFIN-ROE",
    TRUE ~ paste(var1, var2, sep = "-")
  )) %>%
  distinct(year, pair, .keep_all = TRUE)

p_cor_evol <- ggplot(cor_evolution_filtered, aes(x = year, y = cor, color = pair, group = pair)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(title = "Evolución de Correlaciones entre Variables Financieras Clave",
       subtitle = "Cambios en la estructura de relaciones a lo largo del período 2016-2023",
       x = "Año", y = "Correlación", color = "Par de variables") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("evolucion_correlaciones_dual.png", p_cor_evol, width = 10, height = 7, dpi = 300)

# ------------------------------ 16) RESUMEN EN CONSOLA ------------------------------
cat("\n====================================\n")
cat("RESUMEN DEL ANALISIS STATIS DUAL NEUTROSOFICO\n")
cat("====================================\n")

cat("\n--- INTERPRETACION DEL ESPACIO DE COMPROMISO DUAL ---\n")
cat("Dimensión 1 (", round(inercia_pct_dual[1], 2), "%): ")
if (any(coord_vars_dual$Dim1 > 0)) {
  pos_vars <- coord_vars_dual %>% filter(Dim1 > 0) %>% pull(variable) %>% paste(collapse = ", ")
  neg_vars <- coord_vars_dual %>% filter(Dim1 < 0) %>% pull(variable) %>% paste(collapse = ", ")
  cat("Separa (+)", pos_vars, "de (-)", neg_vars, "\n")
}

cat("Dimensión 2 (", round(inercia_pct_dual[2], 2), "%): ")
if (any(coord_vars_dual$Dim2 > 0)) {
  pos_vars <- coord_vars_dual %>% filter(Dim2 > 0) %>% pull(variable) %>% paste(collapse = ", ")
  neg_vars <- coord_vars_dual %>% filter(Dim2 < 0) %>% pull(variable) %>% paste(collapse = ", ")
  cat("Separa (+)", pos_vars, "de (-)", neg_vars, "\n")
}

cat("\n--- TOP VARIABLES POR CONTRIBUCION AL PLANO DUAL ---\n")
print(coord_vars_dual %>% 
        arrange(desc(inercia_contrib)) %>% 
        select(variable, Dim1, Dim2, inercia_contrib, cluster) %>% 
        head(5))

cat("\n--- AÑOS CON MAYOR PESO EN EL COMPROMISO DUAL ---\n")
print(data.frame(year = years, weight = weights_dual) %>% arrange(desc(weight)) %>% head(3))

cat("\n--- AÑOS CON MENOR PESO EN EL COMPROMISO DUAL ---\n")
print(data.frame(year = years, weight = weights_dual) %>% arrange(weight) %>% head(3))

cat("\n====================================\n")
cat("Proceso STATIS Dual Neutrosófico finalizado correctamente.\n")
cat("Archivos generados:\n")
cat("  - coordenadas_variables_compromiso_dual.csv\n")
cat("  - trayectorias_estructura_variables_dual.csv\n")
cat("  - coordenadas_cooperativas_compromiso_dual.csv\n")
cat("  - matriz_rv_dual_neutrosofica.csv\n")
cat("  - pesos_statis_dual_neutrosofico.csv\n")
cat("  - evolucion_correlaciones_dual.png\n")
cat("====================================\n")

# Generar heatmap de incertidumbre para el análisis dual
# (utilizando los mismos datos long del STATIS dual)

I_heat_dual <- long %>%
  group_by(coop, variable) %>%
  summarise(I_mean = mean(I, na.rm = TRUE), .groups = "drop")

p_I_dual <- ggplot(I_heat_dual, aes(x = variable, y = coop, fill = I_mean)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred") +
  labs(title = "Mapa de incertidumbre promedio por cooperativa y variable\n(Modelo STATIS Dual Neutrosófico)",
       x = "Variable", y = "Cooperativa", fill = "I promedio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6))

ggsave("mapa_incertidumbre_dual_neutrosofico.png", p_I_dual, width = 10, height = 12, dpi = 300)


# 15.2.1 Dual compromise map: Dimension 2 vs. Dimension 3

p_vars_dual_dim2_dim3 <- ggplot(
  coord_vars_dual_plot,
  aes(x = Dim2_plot, y = Dim3_plot, color = cluster)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey70"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey70"
  ) +
  
  # Unit correlation circle
  annotate(
    "path",
    x = cos(seq(0, 2 * pi, length.out = 100)),
    y = sin(seq(0, 2 * pi, length.out = 100)),
    color = "grey50",
    linetype = "dotted"
  ) +
  
  # Variable vectors
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = Dim2_plot,
      yend = Dim3_plot
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 0.8,
    color = "darkblue"
  ) +
  
  # Variable points
  geom_point(
    aes(color = cluster),
    size = 4
  ) +
  
  # Variable labels
  ggrepel::geom_text_repel(
    aes(label = variable),
    size = 5,
    fontface = "bold",
    max.overlaps = 100
  ) +
  
  # Titles and axis labels
  labs(
    title = "Neutrosophic Dual STATIS Compromise Map",
    subtitle = sprintf(
      paste0(
        "Consensus structure of relationships among financial variables\n",
        "Inertia Dim2 = %.2f%%, Dim3 = %.2f%%"
      ),
      inercia_pct_dual[2],
      inercia_pct_dual[3]
    ),
    x = paste0(
      "Dimension 2 (",
      round(inercia_pct_dual[2], 2),
      "%)"
    ),
    y = paste0(
      "Dimension 3 (",
      round(inercia_pct_dual[3], 2),
      "%)"
    ),
    color = "Cluster"
  ) +
  
  coord_fixed(ratio = 1, xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# Display the figure
print(p_vars_dual_dim2_dim3)

# Save the figure
ggsave(
  filename = "dual_compromise_map_dim2_dim3.png",
  plot = p_vars_dual_dim2_dim3,
  width = 10,
  height = 8,
  dpi = 300
)

# 15.2.1 Dual compromise map: Dimension 1 vs. Dimension 3

p_vars_dual_dim1_dim3 <- ggplot(
  coord_vars_dual_plot,
  aes(x = Dim1_plot, y = Dim3_plot, color = cluster)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey70"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey70"
  ) +
  
  # Unit correlation circle
  annotate(
    "path",
    x = cos(seq(0, 2 * pi, length.out = 100)),
    y = sin(seq(0, 2 * pi, length.out = 100)),
    color = "grey50",
    linetype = "dotted"
  ) +
  
  # Variable vectors
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = Dim1_plot,
      yend = Dim3_plot
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 0.8,
    color = "darkblue"
  ) +
  
  # Variable points
  geom_point(
    aes(color = cluster),
    size = 4
  ) +
  
  # Variable labels
  ggrepel::geom_text_repel(
    aes(label = variable),
    size = 5,
    fontface = "bold",
    max.overlaps = 100
  ) +
  
  # Titles and axis labels
  labs(
    title = "Neutrosophic Dual STATIS Compromise Map",
    subtitle = sprintf(
      paste0(
        "Consensus structure of relationships among financial variables\n",
        "Inertia Dim1 = %.2f%%, Dim3 = %.2f%%"
      ),
      inercia_pct_dual[1],
      inercia_pct_dual[3]
    ),
    x = paste0(
      "Dimension 1 (",
      round(inercia_pct_dual[1], 2),
      "%)"
    ),
    y = paste0(
      "Dimension 3 (",
      round(inercia_pct_dual[3], 2),
      "%)"
    ),
    color = "Cluster"
  ) +
  
  coord_fixed(ratio = 1, xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# Display the figure
print(p_vars_dual_dim1_dim3)

# Save the figure
ggsave(
  filename = "dual_compromise_map_dim1_dim3.png",
  plot = p_vars_dual_dim1_dim3,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = "Figure_2_Neutrosophic_Dim1_Dim3_CORRECTED.pdf",
  plot = p_vars_dual_dim1_dim3,
  device = cairo_pdf,
  width = 10,
  height = 8,
  units = "in"
)






# ============================================================
# VARIABLE COORDINATES, CONTRIBUTIONS, REPRESENTATION QUALITY,
# AND CLUSTER MEMBERSHIP
# ============================================================

# Names of all retained dimensions
dimension_names <- paste0("Dim", 1:nf)

# Calculate contributions and representation quality
variable_results_table <- coord_vars_dual %>%
  rowwise() %>%
  mutate(
    # Total squared distance using all retained dimensions
    total_squared_distance = sum(
      c_across(all_of(dimension_names))^2,
      na.rm = TRUE
    ),
    
    # Contributions to the first three dimensions
    CTR_Dim1 = 100 * Dim1^2 / lambda_dual[1],
    CTR_Dim2 = 100 * Dim2^2 / lambda_dual[2],
    CTR_Dim3 = 100 * Dim3^2 / lambda_dual[3],
    
    # Representation quality in the Dim1-Dim2 plane
    QLT_2D = 100 * (Dim1^2 + Dim2^2) /
      total_squared_distance,
    
    # Representation quality in the first three dimensions
    QLT_3D = 100 * (Dim1^2 + Dim2^2 + Dim3^2) /
      total_squared_distance
  ) %>%
  ungroup() %>%
  transmute(
    Variable = variable,
    Dim1 = round(Dim1, 3),
    Dim2 = round(Dim2, 3),
    Dim3 = round(Dim3, 3),
    `CTR1 (%)` = round(CTR_Dim1, 2),
    `CTR2 (%)` = round(CTR_Dim2, 2),
    `CTR3 (%)` = round(CTR_Dim3, 2),
    `QLT 2D (%)` = round(QLT_2D, 2),
    `QLT 3D (%)` = round(QLT_3D, 2),
    Cluster = cluster
  )

# Display the complete table
print(variable_results_table)

# Verify the exact k-means cluster assignments
print(
  variable_results_table %>%
    select(Variable, Cluster) %>%
    arrange(Cluster, Variable)
)

# Verify that each dimension's contributions sum to approximately 100%
print(
  variable_results_table %>%
    summarise(
      `Total CTR1` = sum(`CTR1 (%)`),
      `Total CTR2` = sum(`CTR2 (%)`),
      `Total CTR3` = sum(`CTR3 (%)`)
    )
)

# Export the table
write.csv(
  variable_results_table,
  "variable_coordinates_contributions_quality_clusters.csv",
  row.names = FALSE
)









# ============================================================
# CLASSICAL DUAL STATIS ANALYSIS
# Original financial indicators
# ============================================================

# ------------------------------------------------------------
# 1. ANNUAL MATRICES OF ORIGINAL FINANCIAL INDICATORS
# ------------------------------------------------------------

X_years_classical <- lapply(
  years,
  function(current_year) {
    make_year_matrix(
      df_long = long,
      year_value = current_year,
      value_col = "value",
      coops = coops,
      vars_ord = vars
    )
  }
)

names(X_years_classical) <- as.character(years)

# Center and scale each annual matrix separately
X_years_classical_cs <- lapply(
  X_years_classical,
  center_scale_matrix
)

# ------------------------------------------------------------
# 2. CLASSICAL DUAL SCALAR-PRODUCT MATRICES
# ------------------------------------------------------------

W_k_classical <- lapply(
  X_years_classical_cs,
  function(X) {
    t(X) %*% X
  }
)

names(W_k_classical) <- as.character(years)

# ------------------------------------------------------------
# 3. CLASSICAL TEMPORAL RV MATRIX
# ------------------------------------------------------------

K <- length(years)

RV_classical <- matrix(
  0,
  nrow = K,
  ncol = K,
  dimnames = list(
    as.character(years),
    as.character(years)
  )
)

for (i in seq_len(K)) {
  for (j in seq_len(K)) {
    RV_classical[i, j] <- rv_coeff(
      W_k_classical[[i]],
      W_k_classical[[j]]
    )
  }
}

cat("\n====================================\n")
cat("CLASSICAL DUAL STATIS RV MATRIX\n")
cat("====================================\n")
print(round(RV_classical, 6))

# ------------------------------------------------------------
# 4. CLASSICAL COMPROMISE WEIGHTS
# ------------------------------------------------------------

eigen_RV_classical <- eigen(RV_classical)

first_eigenvector_classical <-
  eigen_RV_classical$vectors[, 1]

# Resolve the arbitrary sign of the eigenvector
if (sum(first_eigenvector_classical) < 0) {
  first_eigenvector_classical <-
    -first_eigenvector_classical
}

weights_classical <-
  first_eigenvector_classical /
  sum(first_eigenvector_classical)

names(weights_classical) <- as.character(years)

classical_weights_table <- data.frame(
  Year = years,
  Weight = weights_classical
)

cat("\n====================================\n")
cat("CLASSICAL DUAL STATIS WEIGHTS\n")
cat("====================================\n")
print(
  classical_weights_table %>%
    mutate(Weight = round(Weight, 6))
)

# ------------------------------------------------------------
# 5. CLASSICAL COMPROMISE MATRIX
# ------------------------------------------------------------

W_compromise_classical <- Reduce(
  `+`,
  Map(
    function(Wk, weight) {
      weight * Wk
    },
    W_k_classical,
    weights_classical
  )
)

rownames(W_compromise_classical) <- vars
colnames(W_compromise_classical) <- vars

cat("\n====================================\n")
cat("CLASSICAL DUAL COMPROMISE MATRIX\n")
cat("====================================\n")
print(round(W_compromise_classical, 6))

# ------------------------------------------------------------
# 6. EIGENDECOMPOSITION AND EXPLAINED INERTIA
# ------------------------------------------------------------

eigen_compromise_classical <-
  eigen(W_compromise_classical)

positive_indices_classical <- which(
  eigen_compromise_classical$values > 1e-8
)

nf_classical <- min(
  nf,
  length(positive_indices_classical)
)

lambda_classical <-
  eigen_compromise_classical$values[
    positive_indices_classical
  ][seq_len(nf_classical)]

U_classical <-
  eigen_compromise_classical$vectors[
    ,
    positive_indices_classical,
    drop = FALSE
  ][
    ,
    seq_len(nf_classical),
    drop = FALSE
  ]

rownames(U_classical) <- vars
colnames(U_classical) <-
  paste0("Dim", seq_len(nf_classical))

inertia_classical <- 100 *
  lambda_classical /
  sum(
    eigen_compromise_classical$values[
      eigen_compromise_classical$values > 0
    ]
  )

classical_inertia_table <- data.frame(
  Dimension = paste0(
    "Dim",
    seq_len(nf_classical)
  ),
  Eigenvalue = lambda_classical,
  `Inertia (%)` = inertia_classical,
  `Cumulative inertia (%)` =
    cumsum(inertia_classical),
  check.names = FALSE
)

cat("\n====================================\n")
cat("CLASSICAL COMPROMISE INERTIA\n")
cat("====================================\n")
print(
  classical_inertia_table %>%
    mutate(
      Eigenvalue = round(Eigenvalue, 6),
      `Inertia (%)` =
        round(`Inertia (%)`, 2),
      `Cumulative inertia (%)` =
        round(`Cumulative inertia (%)`, 2)
    )
)

# ------------------------------------------------------------
# 7. CLASSICAL VARIABLE COORDINATES
# ------------------------------------------------------------

coord_vars_classical <- as.data.frame(
  U_classical %*%
    diag(
      sqrt(lambda_classical),
      nrow = nf_classical
    )
)

colnames(coord_vars_classical) <-
  paste0("Dim", seq_len(nf_classical))

coord_vars_classical$variable <- vars

# ------------------------------------------------------------
# 8. CLASSICAL VARIABLE CONTRIBUTIONS
# ------------------------------------------------------------

for (dimension_index in seq_len(nf_classical)) {
  
  dimension_name <-
    paste0("Dim", dimension_index)
  
  contribution_name <-
    paste0("CTR_Dim", dimension_index)
  
  coord_vars_classical[[contribution_name]] <-
    100 *
    coord_vars_classical[[dimension_name]]^2 /
    lambda_classical[dimension_index]
}

# ------------------------------------------------------------
# 9. CLASSICAL VARIABLE CLUSTERING
# ------------------------------------------------------------

set.seed(123)

kmeans_classical <- kmeans(
  coord_vars_classical[
    ,
    c("Dim1", "Dim2")
  ],
  centers = 3,
  nstart = 25
)

coord_vars_classical$cluster <-
  paste0(
    "Cluster ",
    kmeans_classical$cluster
  )

# Reescalado EXCLUSIVAMENTE grafico para los mapas clasicos con circulo unitario.
# Se conserva coord_vars_classical sin ninguna modificacion.
coord_vars_classical_plot <- coord_vars_classical %>%
  mutate(
    Dim1_plot = Dim1 / graph_scale_factor,
    Dim2_plot = Dim2 / graph_scale_factor,
    Dim3_plot = Dim3 / graph_scale_factor
  )

classical_coordinates_table <-
  coord_vars_classical %>%
  transmute(
    Variable = variable,
    Dim1 = round(Dim1, 3),
    Dim2 = round(Dim2, 3),
    Dim3 = round(Dim3, 3),
    Cluster = cluster
  )

cat("\n====================================\n")
cat("CLASSICAL VARIABLE COORDINATES AND CLUSTERS\n")
cat("====================================\n")
print(classical_coordinates_table)

# ------------------------------------------------------------
# 10. CLASSICAL ANNUAL VARIABLE TRAJECTORIES
# ------------------------------------------------------------

partial_classical_coordinates <- lapply(
  seq_along(years),
  function(i) {
    
    inverse_sqrt_lambda <- diag(
      1 / sqrt(lambda_classical),
      nrow = nf_classical
    )
    
    Fk_classical <-
      W_k_classical[[i]] %*%
      U_classical %*%
      inverse_sqrt_lambda
    
    output <- as.data.frame(Fk_classical)
    
    colnames(output) <-
      paste0(
        "Dim",
        seq_len(nf_classical)
      )
    
    output$variable <- vars
    output$year <- years[i]
    output$weight <- weights_classical[i]
    
    output
  }
)

partial_classical_df <-
  bind_rows(partial_classical_coordinates)

# ------------------------------------------------------------
# 11. EXPORT CLASSICAL RESULTS
# ------------------------------------------------------------

write.csv(
  as.data.frame(RV_classical) %>%
    mutate(Year = rownames(RV_classical)) %>%
    relocate(Year),
  "classical_dual_rv_matrix.csv",
  row.names = FALSE
)

write.csv(
  classical_weights_table,
  "classical_dual_weights.csv",
  row.names = FALSE
)

write.csv(
  classical_inertia_table,
  "classical_dual_inertia.csv",
  row.names = FALSE
)

write.csv(
  coord_vars_classical,
  "classical_dual_variable_coordinates.csv",
  row.names = FALSE
)

write.csv(
  partial_classical_df,
  "classical_dual_variable_trajectories.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(W_compromise_classical) %>%
    mutate(Variable = rownames(W_compromise_classical)) %>%
    relocate(Variable),
  "classical_dual_compromise_matrix.csv",
  row.names = FALSE
)









# ============================================================
# CLASSICAL DUAL STATIS COMPROMISE MAP: DIM1 VS. DIM2
# ============================================================

p_classical_dim1_dim2 <- ggplot(
  coord_vars_classical_plot,
  aes(
    x = Dim1_plot,
    y = Dim2_plot,
    color = cluster
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey70"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey70"
  ) +
  annotate(
    "path",
    x = cos(seq(0, 2 * pi, length.out = 100)),
    y = sin(seq(0, 2 * pi, length.out = 100)),
    color = "grey50",
    linetype = "dotted"
  ) +
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = Dim1_plot,
      yend = Dim2_plot
    ),
    arrow = arrow(
      length = unit(0.2, "cm")
    ),
    linewidth = 0.8,
    color = "darkblue"
  ) +
  geom_point(
    size = 4
  ) +
  ggrepel::geom_text_repel(
    aes(label = variable),
    size = 5,
    fontface = "bold",
    max.overlaps = 100
  ) +
  labs(
    title = "Classical Dual STATIS Compromise Map",
    subtitle = sprintf(
      paste0(
        "Consensus structure of relationships among original financial indicators\n",
        "Inertia Dim1 = %.2f%%, Dim2 = %.2f%%"
      ),
      inertia_classical[1],
      inertia_classical[2]
    ),
    x = paste0(
      "Dimension 1 (",
      round(inertia_classical[1], 2),
      "%)"
    ),
    y = paste0(
      "Dimension 2 (",
      round(inertia_classical[2], 2),
      "%)"
    ),
    color = "Cluster"
  ) +
  coord_fixed(ratio = 1, xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

print(p_classical_dim1_dim2)

ggsave(
  filename =
    "classical_dual_compromise_map_dim1_dim2.png",
  plot = p_classical_dim1_dim2,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = "Figure_5_Classical_Dim1_Dim2_CORRECTED.pdf",
  plot = p_classical_dim1_dim2,
  device = cairo_pdf,
  width = 10,
  height = 8,
  units = "in"
)

# ============================================================
# CLASSICAL DUAL STATIS COMPROMISE MAP: DIM1 VS. DIM3
# ============================================================

p_classical_dim1_dim3 <- ggplot(
  coord_vars_classical_plot,
  aes(
    x = Dim1_plot,
    y = Dim3_plot,
    color = cluster
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey70"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey70"
  ) +
  annotate(
    "path",
    x = cos(seq(0, 2 * pi, length.out = 100)),
    y = sin(seq(0, 2 * pi, length.out = 100)),
    color = "grey50",
    linetype = "dotted"
  ) +
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = Dim1_plot,
      yend = Dim3_plot
    ),
    arrow = arrow(
      length = unit(0.2, "cm")
    ),
    linewidth = 0.8,
    color = "darkblue"
  ) +
  geom_point(
    size = 4
  ) +
  ggrepel::geom_text_repel(
    aes(label = variable),
    size = 5,
    fontface = "bold",
    max.overlaps = 100
  ) +
  labs(
    title = "Classical Dual STATIS Compromise Map",
    subtitle = sprintf(
      paste0(
        "Consensus structure of relationships among original financial indicators\n",
        "Inertia Dim1 = %.2f%%, Dim3 = %.2f%%"
      ),
      inertia_classical[1],
      inertia_classical[3]
    ),
    x = paste0(
      "Dimension 1 (",
      round(inertia_classical[1], 2),
      "%)"
    ),
    y = paste0(
      "Dimension 3 (",
      round(inertia_classical[3], 2),
      "%)"
    ),
    color = "Cluster"
  ) +
  coord_fixed(ratio = 1, xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

print(p_classical_dim1_dim3)

ggsave(
  filename =
    "classical_dual_compromise_map_dim1_dim3.png",
  plot = p_classical_dim1_dim3,
  width = 10,
  height = 8,
  dpi = 300
)











# ============================================================
# SENSITIVITY ANALYSIS
# Neutrosophic-score Dual STATIS
# ============================================================

# This block must be executed after the main script.
# Required objects:
# long, years, coops, vars,
# make_year_matrix(), center_scale_matrix(), rv_coeff()

# ------------------------------------------------------------
# 1. SENSITIVITY SCENARIOS
# ------------------------------------------------------------

sensitivity_scenarios <- list(
  
  "Principal specification" = list(
    wT = 0.50,
    wI = 0.30,
    wF = 0.20,
    indeterminacy_scheme = "variable_specific"
  ),
  
  "Equal score weights" = list(
    wT = 1 / 3,
    wI = 1 / 3,
    wF = 1 / 3,
    indeterminacy_scheme = "variable_specific"
  ),
  
  "Greater truth weight" = list(
    wT = 0.60,
    wI = 0.25,
    wF = 0.15,
    indeterminacy_scheme = "variable_specific"
  ),
  
  "Greater determinacy weight" = list(
    wT = 0.40,
    wI = 0.40,
    wF = 0.20,
    indeterminacy_scheme = "variable_specific"
  ),
  
  "Equal indeterminacy weights" = list(
    wT = 0.50,
    wI = 0.30,
    wF = 0.20,
    indeterminacy_scheme = "equal"
  )
)

# ------------------------------------------------------------
# 2. FUNCTION TO RECALCULATE INDETERMINACY
# ------------------------------------------------------------

calculate_scenario_indeterminacy <- function(
    data,
    indeterminacy_scheme = "variable_specific"
) {
  
  if (indeterminacy_scheme == "equal") {
    
    I_scenario <- (
      0.50 * data$vol_local +
        0.50 * data$dev_sector
    )
    
  } else if (indeterminacy_scheme == "variable_specific") {
    
    I_scenario <- dplyr::case_when(
      
      data$variable %in% c("ROA", "ROE") ~
        0.70 * data$vol_local +
        0.30 * data$dev_sector,
      
      data$variable == "MORAMP" ~
        0.60 * data$vol_local +
        0.40 * data$dev_sector,
      
      data$variable == "LIQCOR" ~
        0.50 * data$vol_local +
        0.50 * data$dev_sector,
      
      data$variable %in% c("DEP", "SCC") ~
        0.65 * data$vol_local +
        0.35 * data$dev_sector,
      
      data$variable == "INTFIN" ~
        0.60 * data$vol_local +
        0.40 * data$dev_sector,
      
      TRUE ~
        0.70 * data$vol_local +
        0.30 * data$dev_sector
    )
    
  } else {
    
    stop(
      paste(
        "Unknown indeterminacy scheme:",
        indeterminacy_scheme
      )
    )
  }
  
  pmin(pmax(I_scenario, 0), 1)
}

# ------------------------------------------------------------
# 3. FUNCTION TO RUN DUAL STATIS FOR ONE SCENARIO
# ------------------------------------------------------------

run_sensitivity_scenario <- function(
    data,
    wT,
    wI,
    wF,
    indeterminacy_scheme,
    n_dimensions = 3
) {
  
  if (abs(wT + wI + wF - 1) > 1e-10) {
    stop("The score weights must sum to one.")
  }
  
  scenario_data <- data
  
  scenario_data$I_scenario <-
    calculate_scenario_indeterminacy(
      data = scenario_data,
      indeterminacy_scheme = indeterminacy_scheme
    )
  
  scenario_data$F_scenario <-
    pmin(
      pmax(
        1 - scenario_data$T,
        0
      ),
      1
    )
  
  scenario_data$score_scenario <-
    wT * scenario_data$T +
    wI * (1 - scenario_data$I_scenario) +
    wF * (1 - scenario_data$F_scenario)
  
  # Annual matrices
  X_years_scenario <- lapply(
    years,
    function(current_year) {
      
      make_year_matrix(
        df_long = scenario_data,
        year_value = current_year,
        value_col = "score_scenario",
        coops = coops,
        vars_ord = vars
      )
    }
  )
  
  names(X_years_scenario) <- as.character(years)
  
  # Center and scale annual matrices
  X_years_scenario_cs <- lapply(
    X_years_scenario,
    center_scale_matrix
  )
  
  # Dual scalar-product matrices
  W_years_scenario <- lapply(
    X_years_scenario_cs,
    function(X) {
      t(X) %*% X
    }
  )
  
  names(W_years_scenario) <- as.character(years)
  
  # Temporal RV matrix
  K_scenario <- length(years)
  
  RV_scenario <- matrix(
    0,
    nrow = K_scenario,
    ncol = K_scenario,
    dimnames = list(
      as.character(years),
      as.character(years)
    )
  )
  
  for (i in seq_len(K_scenario)) {
    for (j in seq_len(K_scenario)) {
      
      RV_scenario[i, j] <- rv_coeff(
        W_years_scenario[[i]],
        W_years_scenario[[j]]
      )
    }
  }
  
  # Compromise weights
  eigen_RV_scenario <- eigen(
    RV_scenario,
    symmetric = TRUE
  )
  
  first_eigenvector <-
    eigen_RV_scenario$vectors[, 1]
  
  # Resolve the arbitrary eigenvector sign
  if (sum(first_eigenvector) < 0) {
    first_eigenvector <- -first_eigenvector
  }
  
  weights_scenario <-
    first_eigenvector /
    sum(first_eigenvector)
  
  names(weights_scenario) <- as.character(years)
  
  # Compromise matrix
  W_compromise_scenario <- Reduce(
    `+`,
    Map(
      function(Wk, weight) {
        weight * Wk
      },
      W_years_scenario,
      weights_scenario
    )
  )
  
  rownames(W_compromise_scenario) <- vars
  colnames(W_compromise_scenario) <- vars
  
  # Eigendecomposition
  eigen_compromise_scenario <- eigen(
    W_compromise_scenario,
    symmetric = TRUE
  )
  
  positive_indices <- which(
    eigen_compromise_scenario$values > 1e-8
  )
  
  retained_dimensions <- min(
    n_dimensions,
    length(positive_indices)
  )
  
  lambda_scenario <-
    eigen_compromise_scenario$values[
      positive_indices
    ][seq_len(retained_dimensions)]
  
  U_scenario <-
    eigen_compromise_scenario$vectors[
      ,
      positive_indices,
      drop = FALSE
    ][
      ,
      seq_len(retained_dimensions),
      drop = FALSE
    ]
  
  rownames(U_scenario) <- vars
  
  variable_coordinates <-
    U_scenario %*%
    diag(
      sqrt(lambda_scenario),
      nrow = retained_dimensions
    )
  
  rownames(variable_coordinates) <- vars
  
  colnames(variable_coordinates) <-
    paste0(
      "Dim",
      seq_len(retained_dimensions)
    )
  
  inertia_scenario <-
    100 *
    lambda_scenario /
    sum(
      eigen_compromise_scenario$values[
        eigen_compromise_scenario$values > 1e-8
      ]
    )
  
  list(
    data = scenario_data,
    X_years = X_years_scenario,
    W_years = W_years_scenario,
    RV = RV_scenario,
    weights = weights_scenario,
    compromise = W_compromise_scenario,
    eigenvalues = lambda_scenario,
    inertia = inertia_scenario,
    coordinates = variable_coordinates
  )
}

# ------------------------------------------------------------
# 4. PROCRUSTES CORRELATION
# ------------------------------------------------------------

calculate_procrustes_correlation <- function(
    reference_coordinates,
    comparison_coordinates
) {
  
  X <- as.matrix(reference_coordinates)
  Y <- as.matrix(comparison_coordinates)
  
  common_dimensions <- min(
    ncol(X),
    ncol(Y)
  )
  
  X <- X[, seq_len(common_dimensions), drop = FALSE]
  Y <- Y[, seq_len(common_dimensions), drop = FALSE]
  
  # Center configurations
  X <- scale(
    X,
    center = TRUE,
    scale = FALSE
  )
  
  Y <- scale(
    Y,
    center = TRUE,
    scale = FALSE
  )
  
  norm_X <- sqrt(sum(X^2))
  norm_Y <- sqrt(sum(Y^2))
  
  if (
    norm_X == 0 ||
    norm_Y == 0
  ) {
    return(NA_real_)
  }
  
  X <- X / norm_X
  Y <- Y / norm_Y
  
  singular_values <- svd(
    t(X) %*% Y
  )$d
  
  # Maximum correlation after optimal orthogonal rotation
  as.numeric(
    sum(singular_values)
  )
}

# ------------------------------------------------------------
# 5. RUN ALL SCENARIOS
# ------------------------------------------------------------

sensitivity_results <- lapply(
  sensitivity_scenarios,
  function(specification) {
    
    run_sensitivity_scenario(
      data = long,
      wT = specification$wT,
      wI = specification$wI,
      wF = specification$wF,
      indeterminacy_scheme =
        specification$indeterminacy_scheme,
      n_dimensions = 3
    )
  }
)

principal_result <- sensitivity_results$`Principal specification`

# ------------------------------------------------------------
# 6. CREATE THE SENSITIVITY SUMMARY TABLE
# ------------------------------------------------------------

principal_RV_values <-
  principal_result$RV[
    upper.tri(principal_result$RV)
  ]

sensitivity_summary <- dplyr::bind_rows(
  lapply(
    names(sensitivity_results),
    function(scenario_name) {
      
      result <-
        sensitivity_results[[scenario_name]]
      
      specification <-
        sensitivity_scenarios[[scenario_name]]
      
      scenario_RV_values <-
        result$RV[
          upper.tri(result$RV)
        ]
      
      RV_correlation <- cor(
        principal_RV_values,
        scenario_RV_values,
        method = "pearson"
      )
      
      weight_rank_correlation <- cor(
        principal_result$weights,
        result$weights,
        method = "spearman"
      )
      
      procrustes_correlation <-
        calculate_procrustes_correlation(
          reference_coordinates =
            principal_result$coordinates,
          comparison_coordinates =
            result$coordinates
        )
      
      cumulative_inertia <-
        sum(
          result$inertia[
            seq_len(
              min(
                3,
                length(result$inertia)
              )
            )
          ]
        )
      
      data.frame(
        Scenario = scenario_name,
        
        `Score weights (T, I, F)` =
          sprintf(
            "%.3f, %.3f, %.3f",
            specification$wT,
            specification$wI,
            specification$wF
          ),
        
        `Indeterminacy weights` =
          ifelse(
            specification$indeterminacy_scheme ==
              "variable_specific",
            "Variable-specific",
            "0.50 temporal / 0.50 sectoral"
          ),
        
        `RV-matrix correlation` =
          RV_correlation,
        
        `Weight rank correlation` =
          weight_rank_correlation,
        
        `Procrustes correlation` =
          procrustes_correlation,
        
        `Cumulative inertia Dim1-Dim3 (%)` =
          cumulative_inertia,
        
        check.names = FALSE
      )
    }
  )
)

sensitivity_summary_rounded <-
  sensitivity_summary %>%
  mutate(
    `RV-matrix correlation` =
      round(
        `RV-matrix correlation`,
        4
      ),
    
    `Weight rank correlation` =
      round(
        `Weight rank correlation`,
        4
      ),
    
    `Procrustes correlation` =
      round(
        `Procrustes correlation`,
        4
      ),
    
    `Cumulative inertia Dim1-Dim3 (%)` =
      round(
        `Cumulative inertia Dim1-Dim3 (%)`,
        2
      )
  )

cat("\n")
cat("====================================================\n")
cat("SENSITIVITY ANALYSIS SUMMARY\n")
cat("====================================================\n")

print(
  sensitivity_summary_rounded,
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. EXPORT ANNUAL WEIGHTS FOR ALL SCENARIOS
# ------------------------------------------------------------

sensitivity_weights <- dplyr::bind_rows(
  lapply(
    names(sensitivity_results),
    function(scenario_name) {
      
      data.frame(
        Scenario = scenario_name,
        Year = years,
        Weight = sensitivity_results[[scenario_name]]$weights
      )
    }
  )
)

# Display annual weights
print(sensitivity_weights)


# ------------------------------------------------------------
# 8. EXPORT INERTIA FOR ALL SCENARIOS
# ------------------------------------------------------------

sensitivity_inertia <- dplyr::bind_rows(
  lapply(
    names(sensitivity_results),
    function(scenario_name) {
      
      result <- sensitivity_results[[scenario_name]]
      
      data.frame(
        Scenario = scenario_name,
        Dimension = paste0(
          "Dim",
          seq_along(result$inertia)
        ),
        Eigenvalue = result$eigenvalues,
        `Inertia (%)` = result$inertia,
        `Cumulative inertia (%)` = cumsum(result$inertia),
        check.names = FALSE
      )
    }
  )
)

# Display inertia results
print(sensitivity_inertia)


# ------------------------------------------------------------
# 9. EXPORT VARIABLE COORDINATES
# ------------------------------------------------------------

sensitivity_coordinates <- dplyr::bind_rows(
  lapply(
    names(sensitivity_results),
    function(scenario_name) {
      
      coordinate_data <- as.data.frame(
        sensitivity_results[[scenario_name]]$coordinates
      )
      
      coordinate_data$Variable <- rownames(coordinate_data)
      coordinate_data$Scenario <- scenario_name
      
      coordinate_data %>%
        dplyr::relocate(
          Scenario,
          Variable
        )
    }
  )
)

# Display variable coordinates
print(sensitivity_coordinates)


# ------------------------------------------------------------
# 10. EXPORT RV MATRICES IN LONG FORMAT
# ------------------------------------------------------------

sensitivity_RV_matrices <- dplyr::bind_rows(
  lapply(
    names(sensitivity_results),
    function(scenario_name) {
      
      RV_data <- as.data.frame(
        as.table(
          sensitivity_results[[scenario_name]]$RV
        )
      )
      
      names(RV_data) <- c(
        "Year_1",
        "Year_2",
        "RV"
      )
      
      RV_data$Scenario <- scenario_name
      
      RV_data %>%
        dplyr::relocate(Scenario)
    }
  )
)

# Display the RV matrices in long format
print(sensitivity_RV_matrices)

# ------------------------------------------------------------
# 11. SAVE RESULTS
# ------------------------------------------------------------

write.csv(
  sensitivity_summary_rounded,
  "sensitivity_analysis_summary.csv",
  row.names = FALSE
)

write.csv(
  sensitivity_weights,
  "sensitivity_annual_weights.csv",
  row.names = FALSE
)

write.csv(
  sensitivity_inertia,
  "sensitivity_inertia.csv",
  row.names = FALSE
)

write.csv(
  sensitivity_coordinates,
  "sensitivity_variable_coordinates.csv",
  row.names = FALSE
)

write.csv(
  sensitivity_RV_matrices,
  "sensitivity_RV_matrices.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. VERIFY THAT THE PRINCIPAL SCENARIO MATCHES THE MAIN SCRIPT
# ------------------------------------------------------------

if (exists("RV_dual")) {
  
  cat("\nMaximum absolute difference between the principal")
  cat(" sensitivity RV matrix and the original RV matrix:\n")
  
  print(
    max(
      abs(
        principal_result$RV -
          RV_dual
      )
    )
  )
}

if (exists("weights_dual")) {
  
  cat("\nMaximum absolute difference between the principal")
  cat(" sensitivity weights and the original weights:\n")
  
  print(
    max(
      abs(
        principal_result$weights -
          weights_dual
      )
    )
  )
}

cat("\n")
cat("Sensitivity analysis completed successfully.\n")
cat("Generated files:\n")
cat("  - sensitivity_analysis_summary.csv\n")
cat("  - sensitivity_annual_weights.csv\n")
cat("  - sensitivity_inertia.csv\n")
cat("  - sensitivity_variable_coordinates.csv\n")
cat("  - sensitivity_RV_matrices.csv\n")





























