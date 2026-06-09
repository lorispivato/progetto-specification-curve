
# ==============================================================================
# 1. CARICAMENTO DEL DATASET
# ==============================================================================
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(lme4)
library(lmerTest)

# Caricamento dei dati dal file CSV con l'encodin latin1 per le lettere accentate
df = read.csv(file="data/All_Stimuli_RTs.csv", header=TRUE, sep=",", encoding = "latin1")

# Correzione delle incongruenze nelle etichette della colonna CATEGORY rimuovendo gli spazi iniziali e finali
df$CATEGORY = stringr::str_squish(df$CATEGORY)

# Rinomina alcune etichette degli errori per renderle più chiare e interpretabili
df = df %>%
  mutate(TYPE_OF_ERROR = ifelse(TYPE_OF_ERROR == "", "NO ERROR", TYPE_OF_ERROR))
df = df %>%
  mutate(TYPE_OF_ERROR = ifelse(TYPE_OF_ERROR == "NONE", "NO RESPONSE", TYPE_OF_ERROR))

# Conversione delle variabili categoriali in fattori
df$ID = as.factor(as.character(df$ID))
df$ITEM = as.factor(df$ITEM)
df$CONDITION = as.factor(as.character(df$CONDITION))
df$BLOCK_4 = as.factor(df$BLOCK_4)
df$CATEGORY = as.factor(df$CATEGORY)
df$FREQUENCY_R = as.factor(df$FREQUENCY_R)
df$NEW_PIC = as.factor(df$NEW_PIC)
df$CATEGORY_RESPONSE = as.factor(df$CATEGORY_RESPONSE)
df$TYPE_OF_ERROR = as.factor(df$TYPE_OF_ERROR)

# Creazione del dataframe per codificare l'allineamento basato sul CATEGORY_RESPONSE e CONDITION
mat = data.frame(
  CATEGORY_RESPONSE = c(1, 0, 1, 0),
  CONDITION = c('Basic', 'Basic', 'Category', 'Category'),
  ALIGNMENT = c('non-aligned', 'aligned', 'aligned', 'non-aligned')
)
mat$CATEGORY_RESPONSE = as.factor(mat$CATEGORY_RESPONSE)
mat$CONDITION = as.factor(mat$CONDITION)
mat$ALIGNMENT = as.factor(mat$ALIGNMENT)
df1 = dplyr::inner_join(df, mat)

# Caricamento del Lexique 3.83
lexique = readRDS("data/Lexique383.rds")

# Dal Lexique, si selezionano solo i nomi
lexique_nouns_filtered = lexique %>%
  filter(cgram == "NOM") %>%
  select(ortho, freqfilms2)

# Si effettua un left join tra df e lexique_nouns
df1 = df1 %>%
  left_join(lexique_nouns_filtered, by = c("WORD_FRENCH" = "ortho"))

# Cancellazione delle righe associate con "pile" che hanno frequenza lessicale = 0
df1 = df1 %>%
  filter(!(WORD_FRENCH == "pile" & freqfilms2 == 0))

# Controllo delle righe che non contengono matches
unmatched_rows = df1 %>%
  filter(is.na(freqfilms2))

# Gestione dei casi particolari con frequenza pari a 0
df1$freqfilms2 = ifelse(df1$WORD_FRENCH %in% c("cor français", "thermos"), 0, df1$freqfilms2)

# Gestione dei casi particolari con orthografy differente
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

# Rimozione di oggetti non utilizzati
remove(df, lexique, lexique_nouns_filtered, mat, unmatched_rows, freq_oeil, freq_oeuf, freq_planche)



# ==============================================================================
# 2. PREPARAZIONE DELLE VARIABILI E PULIZIA DATASET
# ==============================================================================
df1$REACTION_TIME     <- as.numeric(df1$REACTION_TIME)
df1$REACTION_TIME_LOG <- log(df1$REACTION_TIME)

# Rimozione delle colonne superflue
df1$AUDIO_FILE         <- NULL
df1$BLOCK_4            <- NULL
df1$FREQUENCY_R        <- NULL
df1$LIST               <- NULL
df1$PICTURE            <- NULL
df1$WORD_FRENCH        <- NULL
df1$NO_RESPONSE        <- NULL
df1$SEM1               <- NULL
df1$RESPONSE           <- NULL
df1$ERROR              <- NULL
df1$ALIGNMENT          <- relevel(df1$ALIGNMENT, ref = "non-aligned")

# Creazione della tabella incrociata
tabella_incrociata <- table(df1$CATEGORY, df1$CATEGORY_RESPONSE)
print(tabella_incrociata)

# Creazione della tabella delle proporzioni
proporzioni <- tabella_incrociata[, 2] / 360
df_plot <- data.frame(
  Category = names(proporzioni),
  Proporzione = as.numeric(proporzioni)
)

# Grafico della tabella delle proporzioni
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

# Creazione della tabella per i multiversi
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
# 5. FUNZIONE DI ESTRAZIONE AGGIORNATA (SOLO BETA PARAMETRICI)
# ==============================================================================
# Creazione della funzione che estrae le informazioni di interesse dai modelli
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
# 6. ESECUZIONE DEI MODELLI (384 CONFIGURAZIONI)
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
# 7. INIZIALIZZAZIONE DELLE SPECIFICATION CURVES SUI BETA
# ==============================================================================

# Iniezione delle colonne fittizie necessarie a specr per evitare glitch nel rendering
df_beta_all <- d_spec %>% mutate(x = "x", y = "y")

# Configurazione forzata del fattore removed_category
df_beta_all$removed_category <- factor(df_beta_all$removed_category, 
                                       levels = c("None", "instruments de musique", "fruits", "professions", "nature", "accessoires"))



# ==============================================================================
# 7. DIVISIONE IN BASE ALLA SCALA DEI COEFFICIENTI BETA
# ==============================================================================

# Sotto-dataset A: SCALA DEI MILLISECONDI REALI (Modelli Lineari m2 e Gamma m4)
beta_scala_ms  <- df_beta_all %>% filter(grepl("m2_", model_id) | grepl("m4_", model_id))

# Sotto-dataset B: SCALA TRASFORMATA (Modelli Log-Linear m1 e InvGauss Log m3)
beta_scala_log <- df_beta_all %>% filter(grepl("m1_", model_id) | grepl("m3_", model_id))



# ==============================================================================
# SOTTO-DATASET A: SCALA MILLISECONDI (Modelli m2 & m4)
# ==============================================================================

# --- MULTIVERSO ESTESO CON LEAVE-ONE-OUT (N = 192) ---

# 4. Effetto Principale ORDER (N=192)
p_ms_384_ord <- plot_specs(beta_scala_ms %>% filter(Interaction == "ORDER"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala MS): Effetto Principale ORDER", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("figures/BETA_MS_384_Order_Pure.png", plot = p_ms_384_ord, width = 4000, height = 4000/1.6, units = "px")

# 5. Interazione ALIGNMENT:CONDITION (N=192)
p_ms_384_disc <- plot_specs(beta_scala_ms %>% filter(Interaction == "ALIGNMENT:CONDITION"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala MS): Interazione ALIGNMENT:CONDITION", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("figures/BETA_MS_384_Alignment_Condition.png", plot = p_ms_384_disc, width = 4000, height = 4000/1.6, units = "px")

# 6. Interazione ORDER:ALIGNMENT (N=192)
p_ms_384_cont <- plot_specs(beta_scala_ms %>% filter(Interaction == "ORDER:ALIGNMENT"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala MS): Interazione ORDER:ALIGNMENT", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("figures/BETA_MS_384_Order_Alignment.png", plot = p_ms_384_cont, width = 4000, height = 4000/1.6, units = "px")



# ==============================================================================
# SOTTO-DATASET B: SCALA TRASFORMATA (Modelli m1 & m3)
# ==============================================================================

# --- MULTIVERSO ESTESO CON LEAVE-ONE-OUT (N = 192) ---

# 10. Effetto Principale ORDER (N=192)
p_log_384_ord <- plot_specs(beta_scala_log %>% filter(Interaction == "ORDER"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala Trasformata): Effetto Principale ORDER", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("figures/BETA_LOG_384_Order_Pure.png", plot = p_log_384_ord, width = 4000, height = 4000/1.6, units = "px")

# 11. Interazione ALIGNMENT:CONDITION (N=192)
p_log_384_disc <- plot_specs(beta_scala_log %>% filter(Interaction == "ALIGNMENT:CONDITION"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala Trasformata): Interazione ALIGNMENT:CONDITION", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("figures/BETA_LOG_384_Alignment_Condition.png", plot = p_log_384_disc, width = 4000, height = 4000/1.6, units = "px")

# 12. Interazione ORDER:ALIGNMENT (N=192)
p_log_384_cont <- plot_specs(beta_scala_log %>% filter(Interaction == "ORDER:ALIGNMENT"), choices = c("modeling", "random_effects", "datasets", "removed_category")) +
  plot_annotation(title = "Beta Curve (N=192, Scala Trasformata): Interazione ORDER:ALIGNMENT", 
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("figures/BETA_LOG_384_Order_Alignment.png", plot = p_log_384_cont, width = 4000, height = 4000/1.6, units = "px")




# ==============================================================================
# 8. SALVATAGGIO DATI
# ==============================================================================

# save.image("data/workspace.Rdata")
load("data/workspace.Rdata")



# ==============================================================================
# 9. VOLCANO PLOT
# ==============================================================================

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
    
    nome_file = paste0("figures/volcano_", str_replace_all(which_interaction, ":", "_"), "_", 
                       which_scala, ".png")
    ggsave(nome_file, plot = grafico, 
           width = 4000, height = 4000/1.6, units = "px")
  }
}


