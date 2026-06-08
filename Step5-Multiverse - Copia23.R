###########################

#### Step 1: Data loading ####

library(dplyr)
library(stringr)
library(pals)             # For color palettes
# library(ggpubfigs)        # For color palettes
library(RColorBrewer)     # For color palettes
library(ggplot2)          # For creating visualizations
library(ggrepel)          # For geom_text_repel
library(tidyr)            # For data manipulation
library(sjPlot)           # For additional plotting functions
library(PupillometryR)    # For geom_flat_violin
library(patchwork)        # For combining multiple plots
library(lme4)             # For fitting linear mixed-effects models
library(lmerTest)

# Load data from CSV file with Latin1 encoding for accented letters
df = read.csv(file="data/All_Stimuli_RTs.csv", header=TRUE, sep=",", encoding = "latin1")

# Fix label inconsistency in CATEGORY column by removing leading and trailing spaces
df$CATEGORY = stringr::str_squish(df$CATEGORY)

# Rename some labels for errors to make them clearer
df = df %>%
  mutate(TYPE_OF_ERROR = ifelse(TYPE_OF_ERROR == "", "NO ERROR", TYPE_OF_ERROR))
df = df %>%
  mutate(TYPE_OF_ERROR = ifelse(TYPE_OF_ERROR == "NONE", "NO RESPONSE", TYPE_OF_ERROR))

# Convert variables to factors for categorical analysis
df$ID = as.factor(as.character(df$ID))                   # ID as factor
df$ITEM = as.factor(df$ITEM)                             # ITEM as factor
df$CONDITION = as.factor(as.character(df$CONDITION))     # CONDITION as factor
df$BLOCK_4 = as.factor(df$BLOCK_4)                       # BLOCK_4 as factor
df$CATEGORY = as.factor(df$CATEGORY)                     # CATEGORY as factor
df$FREQUENCY_R = as.factor(df$FREQUENCY_R)               # FREQUENCY_R as factor
df$NEW_PIC = as.factor(df$NEW_PIC)                       # NEW_PIC as factor
df$CATEGORY_RESPONSE = as.factor(df$CATEGORY_RESPONSE)   # CATEGORY_RESPONSE as factor
df$TYPE_OF_ERROR = as.factor(df$TYPE_OF_ERROR)           # TYPE_OF_ERROR as factor

# Create a data frame to encode alignment based on CATEGORY_RESPONSE and CONDITION
mat = data.frame(
  CATEGORY_RESPONSE = c(1, 0, 1, 0),
  CONDITION = c('Basic', 'Basic', 'Category', 'Category'),
  ALIGNMENT = c('non-aligned', 'aligned', 'aligned', 'non-aligned')
)
mat$CATEGORY_RESPONSE = as.factor(mat$CATEGORY_RESPONSE)  # Convert CATEGORY_RESPONSE to factor
mat$CONDITION = as.factor(mat$CONDITION)                  # Convert CONDITION to factor
mat$ALIGNMENT = as.factor(mat$ALIGNMENT)                  # Convert ALIGNMENT to factor

# Merge data frames to include alignment information in df1
df1 = dplyr::inner_join(df, mat)

# Load Lexique 3.83
lexique = readRDS("data/Lexique383.rds")

# From Lexique, select only nouns
lexique_nouns_filtered = lexique %>%
  filter(cgram == "NOM") %>%
  select(ortho, freqfilms2)

# Perform a left join between df and lexique_nouns
df1 = df1 %>%
  left_join(lexique_nouns_filtered, by = c("WORD_FRENCH" = "ortho"))

# Delete the rows associated with "pile" that have lexical frequency = 0
df1 = df1 %>%
  filter(!(WORD_FRENCH == "pile" & freqfilms2 == 0))

# Check for rows where no match was found
unmatched_rows = df1 %>%
  filter(is.na(freqfilms2))

# Handle custom cases where freqfilms2 should be set to 0
df1$freqfilms2 = ifelse(df1$WORD_FRENCH %in% c("cor français", "thermos"), 0, df1$freqfilms2)

# Handle custom cases where orthography is different
freq_planche = lexique_nouns_filtered %>% filter(ortho == "planche") %>% pull(freqfilms2)
freq_oeil = lexique_nouns_filtered %>% filter(ortho == "oeil") %>% pull(freqfilms2)
freq_oeuf = lexique_nouns_filtered %>% filter(ortho == "oeuf") %>% pull(freqfilms2)

df1 = df1 %>%
  mutate(
    freqfilms2 = case_when(
      WORD_FRENCH == "planche à découper" ~ freq_planche,  # Set to value of "planche"
      WORD_FRENCH == "œil" ~ freq_oeil,                    # Set to value of "oeil"
      WORD_FRENCH == "œuf" ~ freq_oeuf,                    # Set to value of "oeuf"
      TRUE ~ freqfilms2  # Keep existing freqfilms2 for all other cases
    )
  )

# Rename freqfilms2 to FREQUENCY_STIM
df1 = df1 %>%
  rename(FREQUENCY_STIM = freqfilms2)

# Add a log frequency
df1$FREQUENCY_STIM_LOG = log1p(df1$FREQUENCY_STIM)

# Remove intermediate objects that I don't need anymore
remove(lexique, lexique_nouns_filtered, mat, unmatched_rows, freq_oeil, freq_oeuf, freq_planche)


# ==============================================================================
# 1. IMPOSTAZIONE PALETTE COLORI PER LE CATEGORIE
# ==============================================================================
ordered_categories = as.character(sort(unique(df1$CATEGORY)))
category_colors    <- setNames(stepped3(length(ordered_categories)), ordered_categories)

custom_scale_fill  <- scale_fill_manual(values = category_colors, drop = FALSE)
custom_scale_color <- scale_color_manual(values = category_colors, drop = FALSE)


# ==============================================================================
# 2. PREPARAZIONE DELLE VARIABILI E PULIZIA DATASET
# ==============================================================================
df1$REACTION_TIME     <- as.numeric(df1$REACTION_TIME)
df1$REACTION_TIME_LOG <- log(df1$REACTION_TIME)

# Rimozione colonne superflue
df1$AUDIO_FILE         <- NULL
df1$BLOCK_4            <- NULL
df1$FREQUENCY_R        <- NULL
df1$LIST               <- NULL
df1$PICTURE            <- NULL
df1$WORD_FRENCH        <- NULL
df1$NO_RESPONSE        <- NULL
df1$SEM1               <- NULL
df1$FREQUENCY_STIM     <- NULL
df1$FREQUENCY_STIM_LOG <- NULL
df1$RESPONSE           <- NULL
df1$ERROR              <- NULL

df1$ALIGNMENT          <- relevel(df1$ALIGNMENT, ref = "non-aligned")

# Crea la tabella incrociata
tabella_incrociata <- table(df1$CATEGORY, df1$CATEGORY_RESPONSE)
print(tabella_incrociata)

proporzioni <- tabella_incrociata[, 2] / 360
df_plot <- data.frame(
  Category = names(proporzioni),
  Proporzione = as.numeric(proporzioni)
)

ggplot(df_plot, aes(x = reorder(Category, -Proporzione), y = Proporzione)) +
  geom_col(width = 0.6, fill = "skyblue", color = "black", linewidth = 0.2) + 
  geom_text(aes(label = round(Proporzione, 3)), vjust = -0.5, fontface = "plain", size = 3) + 
  labs(
    title = "Proporzione di risposte Category",
    x = "Category",
    y = "Proportion of Category response"
  ) +
  theme_minimal(base_size = 13) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 5)),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  ylim(0, max(df_plot$Proporzione) + 0.05)


# ==============================================================================
# 3. MULTIVERSE - DEFINIZIONE DELLE OPERAZIONI SUI DATI (DECISION GRID)
# ==============================================================================
scelte_immagini <- c("All_Pictures", "No_Novel_Pictures")
scelte_errori   <- c("Strict_Errors_Only", "Broad_Errors_Included")
scelte_rimozione_cat <- c("None", "instruments de musique", "fruits", "professions", "nature", "accessoires")

griglia_multiverso <- expand.grid(
  Immagini     = scelte_immagini, 
  Errori       = scelte_errori, 
  RimozioneCat = scelte_rimozione_cat,
  stringsAsFactors = FALSE
)


# ==============================================================================
# 4. CONFIGURAZIONE AMBIENTE PER I MODELLI MISTI
# ==============================================================================
options(scipen = 999)
library(broom)
library(broom.mixed)
library(purrr)
library(specr)


# ==============================================================================
# 6. FUNZIONE DI ESTRAZIONE AGGIORNATA (SOLO BETA PARAMETRICI)
# ==============================================================================
extract_for_sca <- function(model, model_name, scelta_img, scelta_err, scelta_cat) {
  
  family_choice <- case_when(
    grepl("m1_", model_name) ~ "Log-Linear",
    grepl("m2_", model_name) ~ "Linear",
    grepl("m3_", model_name) ~ "InvGauss-Log",
    grepl("m4_", model_name) ~ "Gamma-Identity"
  )
  
  group_choice <- if_else(grepl("_item", model_name), "Item", "Category")
  rs_choice    <- if_else(grepl("_rs", model_name), "RandomSlope_Sì", "RandomSlope_No")
  
  # Estrazione dei Beta fissi dal sommario nativo del modello
  beta_tidy <- tidy(model, effects = "fixed") %>% as.data.frame()
  if (!"p.value" %in% colnames(beta_tidy)) beta_tidy$p.value <- NA
  
  beta_raw <- beta_tidy %>%
    filter(
      term == "ORDER" |
        (grepl("ALIGNMENT", term, ignore.case = TRUE) & grepl("CONDITION", term, ignore.case = TRUE)) |
        (grepl("ORDER", term, ignore.case = TRUE) & grepl("ALIGNMENT", term, ignore.case = TRUE))
    ) %>% 
    mutate(Interaction = case_when(
      term == "ORDER" ~ "ORDER",
      grepl("ORDER", term, ignore.case = TRUE) ~ "ORDER:ALIGNMENT",
      TRUE ~ "ALIGNMENT:CONDITION"
    )) %>%
    group_by(Interaction) %>% slice(1) %>% ungroup()
  
  beta_raw %>%
    select(Interaction, estimate = estimate, std.error = std.error, p.value = p.value) %>%
    mutate(
      metric    = "Beta",
      conf.low  = estimate - (1.96 * std.error), 
      conf.high = estimate + (1.96 * std.error),
      model_id         = model_name, 
      modeling         = family_choice, 
      random_effects   = paste("group", group_choice, rs_choice, sep = "_"), 
      datasets         = paste(scelta_img, scelta_err, sep = " + "), 
      removed_category = scelta_cat 
    )
}


# ==============================================================================
# 5 & 7. AUTOMAZIONE MASSIFICA: ECO-CICLO MODELLI (384 CONFIGURAZIONI)
# ==============================================================================
d_spec_accumulatore <- list()

for (i in 1:nrow(griglia_multiverso)) {
  
  riga_corrente <- griglia_multiverso[i, ]
  
  cat("\n==================================================================")
  cat("\nCella Data-Multiverso", i, "di 24")
  cat("\nImmagini:", riga_corrente$Immagini, "| Errori:", riga_corrente$Errori)
  cat("\nCategoria Rimossa:", riga_corrente$RimozioneCat)
  cat("\n==================================================================\n")
  
  dati <- df1
  
  if (riga_corrente$Immagini == "No_Novel_Pictures") {
    dati <- dati %>% filter(NEW_PIC == 0)
  }
  
  if (riga_corrente$Errori == "Strict_Errors_Only") {
    dati <- dati %>% filter(TYPE_OF_ERROR %in% c("NO ERROR", "SYN", "CAT"))
  } else {
    dati <- dati %>% filter(TYPE_OF_ERROR %in% c("NO ERROR", "SYN", "SPE", "SEM1", "CAT"))
  }
  
  if (riga_corrente$RimozioneCat != "None") {
    dati <- dati %>% filter(CATEGORY != riga_corrente$RimozioneCat)
  }
  
  if (nrow(dati) == 0) next
  
  dati$ID        <- as.factor(as.character(dati$ID))
  dati$ITEM      <- as.factor(as.character(dati$ITEM))
  dati$CATEGORY  <- as.factor(as.character(dati$CATEGORY))
  dati$ALIGNMENT <- as.factor(as.character(dati$ALIGNMENT))
  dati$ALIGNMENT <- relevel(dati$ALIGNMENT, ref = "non-aligned")
  
  # Stime Blocco A (No Random Slope)
  m1_item <- lmer(REACTION_TIME_LOG ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1|ID) + (1|ITEM), data = dati, REML = FALSE)
  m2_item <- lmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1|ID) + (1|ITEM), data = dati, REML = FALSE)
  m3_item <- glmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1|ID) + (1|ITEM), data = dati, family = inverse.gaussian(link = log), nAGQ = 0)
  m4_item <- glmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1|ID) + (1|ITEM), data = dati, family = Gamma(link = identity), nAGQ = 0)
  
  m1_category <- lmer(REACTION_TIME_LOG ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1|ID) + (1|CATEGORY), data = dati, REML = FALSE)
  m2_category <- lmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1|ID) + (1|CATEGORY), data = dati, REML = FALSE)
  m3_category <- glmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1|ID) + (1|CATEGORY), data = dati, family = inverse.gaussian(link = log), nAGQ = 0)
  m4_category <- glmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1|ID) + (1|CATEGORY), data = dati, family = Gamma(link = identity), nAGQ = 0)
  
  # Stime Blocco B (Con Random Slope)
  m1_item_rs <- lmer(REACTION_TIME_LOG ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1+CONDITION|ID) + (1|ITEM), data = dati, REML = FALSE)
  m2_item_rs <- lmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1+CONDITION|ID) + (1|ITEM), data = dati, REML = FALSE)
  m3_item_rs <- glmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1+CONDITION|ID) + (1|ITEM), data = dati, family = inverse.gaussian(link = log), nAGQ = 0)
  m4_item_rs <- glmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1+CONDITION|ID) + (1|ITEM), data = dati, family = Gamma(link = identity), nAGQ = 0)
  
  m1_category_rs <- lmer(REACTION_TIME_LOG ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1+CONDITION|ID) + (1|CATEGORY), data = dati, REML = FALSE)
  m2_category_rs <- lmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1+CONDITION|ID) + (1|CATEGORY), data = dati, REML = FALSE)
  m3_category_rs <- glmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1+CONDITION|ID) + (1|CATEGORY), data = dati, family = inverse.gaussian(link = log), nAGQ = 0)
  m4_category_rs <- glmer(REACTION_TIME ~ (ORDER + ALIGNMENT + CONDITION)^2 + (1+CONDITION|ID) + (1|CATEGORY), data = dati, family = Gamma(link = identity), nAGQ = 0)
  
  modelli_correnti <- list(
    "m1_item" = m1_item, "m2_item" = m2_item, "m3_item" = m3_item, "m4_item" = m4_item,
    "m1_category" = m1_category, "m2_category" = m2_category, "m3_category" = m3_category, "m4_category" = m4_category,
    "m1_item_rs" = m1_item_rs, "m2_item_rs" = m2_item_rs, "m3_item_rs" = m3_item_rs, "m4_item_rs" = m4_item_rs,
    "m1_category_rs" = m1_category_rs, "m2_category_rs" = m2_category_rs, "m3_category_rs" = m3_category_rs, "m4_category_rs" = m4_category_rs
  )
  
  d_spec_cella <- lapply(names(modelli_correnti), function(nome) {
    extract_for_sca(modelli_correnti[[nome]], nome, riga_corrente$Immagini, riga_corrente$Errori, riga_corrente$RimozioneCat)
  }) %>% bind_rows()
  
  d_spec_accumulatore[[i]] <- d_spec_cella
  
  rm(list = c(names(modelli_correnti), "modelli_correnti", "dati"))
  gc()
}

# Unione globale di tutti i coefficienti beta
d_spec <- bind_rows(d_spec_accumulatore)


# ==============================================================================
# 10. GENERAZIONE SPECIFICATION CURVES SUI BETA (SCALATI SEPARATAMENTE)
# ==============================================================================

# Iniezione delle colonne fittizie necessarie a specr per evitare glitch nel rendering
df_beta_all <- d_spec %>% mutate(x = "x", y = "y")

# Configurazione forzata del fattore removed_category
df_beta_all$removed_category <- factor(df_beta_all$removed_category, 
                                       levels = c("None", "instruments de musique", "fruits", "professions", "nature", "accessoires"))

# ------------------------------------------------------------------------------
# DIVISIONE IN BASE ALLA SCALA DEI COEFFICIENTI BETA
# ------------------------------------------------------------------------------
# Sotto-dataset A: SCALA DEI MILLISECONDI REALI (Modelli Lineari m2 e Gamma m4)
beta_scala_ms  <- df_beta_all %>% filter(grepl("m2_", model_id) | grepl("m4_", model_id))

# Sotto-dataset B: SCALA TRASFORMATA (Modelli Log-Linear m1 e InvGauss Log m3)
beta_scala_log <- df_beta_all %>% filter(grepl("m1_", model_id) | grepl("m3_", model_id))

# Creazione dei sotto-dataset per il Multiverso Compatto (N = 32 per ciascuna scala)
mv_64_ms  <- beta_scala_ms %>% filter(removed_category == "None")
mv_64_log <- beta_scala_log %>% filter(removed_category == "None")


# ==============================================================================
# BLOCCO SALVATAGGI 1: SCALA MILLISECONDI (Modelli m2 & m4)
# ==============================================================================

# --- MULTIVERSO ESTESO CON LEAVE-ONE-OUT (N = 192) ---

# 4. Effetto Principale ORDER (N=192)
p_ms_384_ord <- plot_specs(beta_scala_ms %>% filter(Interaction == "ORDER"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala MS): Effetto Principale ORDER", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("BETA_MS_384_Order_Pure.png", plot = p_ms_384_ord, width = 4000, height = 4000/1.6, units = "px")

# 5. Interazione ALIGNMENT:CONDITION (N=192)
p_ms_384_disc <- plot_specs(beta_scala_ms %>% filter(Interaction == "ALIGNMENT:CONDITION"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala MS): Interazione ALIGNMENT:CONDITION", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("BETA_MS_384_Alignment_Condition.png", plot = p_ms_384_disc, width = 4000, height = 4000/1.6, units = "px")

# 6. Interazione ORDER:ALIGNMENT (N=192)
p_ms_384_cont <- plot_specs(beta_scala_ms %>% filter(Interaction == "ORDER:ALIGNMENT"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala MS): Interazione ORDER:ALIGNMENT", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("BETA_MS_384_Order_Alignment.png", plot = p_ms_384_cont, width = 4000, height = 4000/1.6, units = "px")


# ==============================================================================
# BLOCCO SALVATAGGI 2: SCALA TRASFORMATA (Modelli m1 & m3)
# ==============================================================================

# --- MULTIVERSO ESTESO CON LEAVE-ONE-OUT (N = 192) ---

# 10. Effetto Principale ORDER (N=192)
p_log_384_ord <- plot_specs(beta_scala_log %>% filter(Interaction == "ORDER"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala Trasformata): Effetto Principale ORDER", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("BETA_LOG_384_Order_Pure.png", plot = p_log_384_ord, width = 4000, height = 4000/1.6, units = "px")

# 11. Interazione ALIGNMENT:CONDITION (N=192)
p_log_384_disc <- plot_specs(beta_scala_log %>% filter(Interaction == "ALIGNMENT:CONDITION"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala Trasformata): Interazione ALIGNMENT:CONDITION", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("BETA_LOG_384_Alignment_Condition.png", plot = p_log_384_disc, width = 4000, height = 4000/1.6, units = "px")

# 12. Interazione ORDER:ALIGNMENT (N=192)
p_log_384_cont <- plot_specs(beta_scala_log %>% filter(Interaction == "ORDER:ALIGNMENT"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala Trasformata): Interazione ORDER:ALIGNMENT", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("BETA_LOG_384_Order_Alignment.png", plot = p_log_384_cont, width = 4000, height = 4000/1.6, units = "px")


save.image("workspace.Rdata")
load("workspace.Rdata")


# Volcano plot -------------------------------------------------------------

# install.packages("BiocManager")
# BiocManager::install("EnhancedVolcano")

# library(EnhancedVolcano)
# 
# df <- beta_scala_ms %>%
#   filter(Interaction == "ORDER:ALIGNMENT") %>%
#   as.data.frame()
# 
# EnhancedVolcano(
#   df,
#   lab = rownames(df),
#   x = "estimate",
#   y = "p.value"
# )



## Scala dei ms ------------------------------------------------------------

beta_all = beta_scala_ms %>% 
  add_row(beta_scala_log) %>% 
  mutate(scala = c(rep("ms", nrow(beta_scala_ms)),
                   rep("log", nrow(beta_scala_log))))

which_interaction = "ORDER"
which_scala = "log"

beta_all %>% 
  filter(scala == which_scala) %>% 
  filter(Interaction == which_interaction) %>% 
  mutate(significativo = (p.value < 0.05)) %>% 
  mutate(p.value = -log10(p.value)) %>% 
  arrange(estimate) %>% 
  ggplot(aes(x=estimate, y=p.value, pch=removed_category, 
         col=significativo)) + 
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col=2, lwd=1) +
  annotate("text", x = -Inf, y = -log10(0.05), label = "Soglia 0.05", 
           vjust = -0.5, hjust = -0.1) +
  theme_bw() +
  labs(x = "Stima del coefficiente", y = "-log10(p-value)",
       pch = "Categoria rimossa", col = "Significativo",
       title = paste0("Volcano plot per ", which_interaction, ", scala ", which_scala)) +
  scale_color_manual(
    values = c("FALSE" = "#F8766D", "TRUE" = "#00BFC4"),
    labels = c("FALSE" = "No", "TRUE" = "Sì")
  ) +
  guides(
    colour = guide_legend(order = 1),
    shape  = guide_legend(order = 2)
  )


# Salvataggio automatico di tutti i 6 grafici

for(which_scala in c("ms", "log")) {
  for(which_interaction in c("ORDER:ALIGNMENT", "ALIGNMENT:CONDITION", "ORDER")) {
    print(c(which_scala, which_interaction))
    
    grafico = beta_all %>% 
      filter(scala == which_scala) %>% 
      filter(Interaction == which_interaction) %>% 
      mutate(significativo = (p.value < 0.05)) %>% 
      mutate(p.value = -log10(p.value)) %>% 
      arrange(estimate) %>% 
      ggplot(aes(x=estimate, y=p.value, pch=removed_category, 
                 col=significativo)) + 
      geom_point() +
      geom_hline(yintercept = -log10(0.05), col=2, lwd=1) +
      annotate("text", x = -Inf, y = -log10(0.05), label = "Soglia 0.05", 
               vjust = -0.5, hjust = -0.1) +
      theme_bw() +
      labs(x = "Stima del coefficiente", y = "-log10(p-value)",
           pch = "Categoria rimossa", col = "Significativo",
           title = paste0("Volcano plot per ", which_interaction, ", scala ", which_scala)) +
      scale_color_manual(
        values = c("FALSE" = "#F8766D", "TRUE" = "#00BFC4"),
        labels = c("FALSE" = "No", "TRUE" = "Sì")
      ) +
      guides(
        colour = guide_legend(order = 1),
        shape  = guide_legend(order = 2)
      )
    
    nome_file = paste0("volcano_", str_replace_all(which_interaction, ":", "_"), "_", 
                       which_scala, ".png")
    ggsave(nome_file, plot = grafico, 
           width = 4000, height = 4000/1.6, units = "px")
  }
}





