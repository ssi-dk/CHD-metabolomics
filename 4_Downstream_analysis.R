# Statistical Analysis for Metabolomics CHD paper

library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(survival)  
library(stringr)
library(readr)
library(forcats)
library(mice)
library(gt)

set.seed(1234)
# ---------------------------------------------------------------------------
# Paths — edit these before running
# ---------------------------------------------------------------------------

dir_main   <- "path/to/main_dir"
dir_input  <- "path/to/data"
dir_output <- "path/to/output"

# Feature table with metabolomics data
ft <- read.csv(paste0(dir_input, "chd_quant_all_cont_blankrm_wavenorm.csv"))

# Metadata containing sample information, including CHD case/control status
md <- read.csv(paste0(dir_input, "metadata.csv"))

# Metabolite annotations
annotations <- read.xlsx(paste0(dir_input, "merged_iroa_gnps_sirius_curated.xlsx"))

# here we remove features related to an MS contaminant, a molecular family with very 
# early eluting features (0.3-0.6) with a
# GNPS library match to MS_Contaminant_Sodium_Formate_Cluster

ms_contaminants <- c("X105","X110","X113","X126","X128","X129","X137","X204","X222","X229","X522","X539","X544","X551","X584","X585","X586","X588","X589","X590","X592","X593","X600","X603","X605","X608","X610","X613","X617","X624","X627","X642","X645","X652","X653","X662","X664","X667","X668","X670","X671","X673","X677","X679","X681","X699","X81","X86","X99")
ft <- ft[, !colnames(ft) %in% ms_contaminants]

# Re-code case/control assignments to 1/0
md$CHD <- ifelse(md$CASE.CONTROL=="CASE",1,0) 

#Check the tables
md[1:10,1:10]
ft[1:10,1:10]

# Merge metadata and feature table
Data <- merge(md,ft,by="SampleName")

Data[1:10,1:10]

dim(Data)
dim(md)
dim(ft)

annotations <- annotations %>% dplyr::select(feature_id,curated_annotation,csi_name)
colnames(annotations)[1] <- "metID"

#Metabolite IDs start with X in the R data tables for simplicity
annotations$metID <- paste0("X",as.character(annotations$metID))

Data$indikatorvariabel <- substr(Data$indikatorvariabel,1,4)

# Remove outliers i.e. samples that either did not contain EDTA or had a low volume
sum(Data$anticoagulant_outlier)
sum(Data$low_vol)
sum(Data$hemolysis)

p_GD_sample <- ggplot(Data,aes(x=as.factor(anticoagulant_outlier),y=GD_sample,color=as.factor(anticoagulant_outlier)))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(alpha=0.5,width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none")+
  xlab("Anticoagulant outlier")+
  ylab("GD_sample")+
  scale_color_jama()

ggsave(paste0(dir_output, "p_GD_sample.pdf"),p_GD_sample)

dim(Data)

# anticoagulant and volume outliers are removed
Data <- subset(Data,anticoagulant_outlier==0  & low_vol==0) 
dim(Data)

# we want to start from index of the column that is after Batch -- that's where the PCs are going to be added
start <- which(names(Data)=="Batch")+1

Data_pairs <- Data[Data$indikatorvariabel %in% Data$indikatorvariabel[duplicated(Data$indikatorvariabel)],]
dim(Data_pairs)
table(Data_pairs$indikatorvariabel)

out_dir <- paste0(dir_output,"results_metab_clogit")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# Tidy covariates & keep complete matched pairs
Data1 <- Data_pairs %>%
  mutate(
    CHD = as.integer(CHD),
    indikatorvariabel = as.factor(indikatorvariabel),
    
    KOEN_BARN = factor(KOEN_BARN),
    BMI = suppressWarnings(as.numeric(BMI)),
    GESTATIONSALDER_DAGE = suppressWarnings(as.numeric(GESTATIONSALDER_DAGE)),
    
    RYGERSTATUS_chr = as.character(smoking),
    RYGERSTATUS_chr = na_if(RYGERSTATUS_chr, "90"),
    RYGERSTATUS_chr = na_if(RYGERSTATUS_chr, "99"),
    RYGERSTATUS_chr = ifelse(RYGERSTATUS_chr %in% c("0","1","2","3"), RYGERSTATUS_chr, NA_character_),
    RYG_cat = factor(case_when(
      RYGERSTATUS_chr == "0" ~ "No",
      RYGERSTATUS_chr == "1" ~ "Occasional",
      RYGERSTATUS_chr == "2" ~ "Daily <15",
      RYGERSTATUS_chr == "3" ~ "Daily ≥15",
      TRUE ~ NA_character_
    ), levels = c("No","Occasional","Daily <15","Daily ≥15")),
    
    alcohol = factor(alcohol, levels = as.character(seq(0, 5, 0.5)), ordered = TRUE),
    parity  = factor(parity),
    
    CHD_subtype = factor(CHD_type1),
  ) 

# Add an explicit row id to merge back later
Data1 <- Data1 %>%
  mutate(.row_id = dplyr::row_number())

# Subset of variables to feed into mice:
# CHD + indikatorvariabel: used as predictors, NOT imputed
# KOEN_BARN, BMI, GESTATIONSALDER_DAGE, RYG_cat, alcohol, parity: to be imputed
imp_dat <- Data1 %>%
  select(.row_id,
         CHD,
         indikatorvariabel,
         KOEN_BARN,
         BMI,
         GESTATIONSALDER_DAGE,
         RYG_cat,
         alcohol,
         parity)

# Initialize mice (gets default method + predictor matrix)
ini  <- mice(imp_dat, maxit = 0, print = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

# NO imputation for CHD or indikatorvariabel
meth[c("CHD", "indikatorvariabel", ".row_id")] <- ""
pred[ , c("CHD", "indikatorvariabel", ".row_id")] <- 0  # they don't use themselves as predictors

# Choose methods per variable type
meth["KOEN_BARN"]            <- "logreg"   # binary factor
meth["BMI"]                  <- "pmm"      # continuous
meth["GESTATIONSALDER_DAGE"] <- "pmm"      # continuous
meth["RYG_cat"]              <- "polyreg"  # unordered categorical
meth["alcohol"]              <- "polr"     # ordered factor
meth["parity"]               <- "polyreg"  # unordered categorical

imp <- mice(
  imp_dat,
  method          = meth,
  predictorMatrix = pred,
  m               = 1,      #single imputation
  maxit           = 20,
  print           = FALSE
)

# Extract completed data
imp_complete <- complete(imp, 1)

# Join imputed covariates back into Data1
Data1 <- Data1 %>%
  select(-KOEN_BARN, -BMI, -GESTATIONSALDER_DAGE,
         -RYG_cat, -alcohol, -parity) %>%
  left_join(
    imp_complete %>%
      select(.row_id,
             KOEN_BARN,
             BMI,
             GESTATIONSALDER_DAGE,
             RYG_cat,
             alcohol,
             parity),
    by = ".row_id"
  )

colSums(is.na(Data1[, c("KOEN_BARN","BMI","GESTATIONSALDER_DAGE",
                        "RYG_cat","alcohol","parity")]))


# Check matched pairs (should be size 2 per stratum)
pairs_ok <- Data1 %>%
  count(indikatorvariabel, name = "n") %>%
  summarise(all_n2 = all(n == 2)) %>% pull(all_n2)
if (!isTRUE(pairs_ok)) {
  warning("Some matched strata are not size 2 after cleaning; they will be ignored by clogit automatically.")
}

# All metabolites are named X\d+
met_names <- met_names <- names(Data1)[grepl("^X[0-9]+$", names(Data1))]

# Coerce metabolites to numeric
Data1[met_names] <- lapply(Data1[met_names], function(x) suppressWarnings(as.numeric(x)))

# Keep only informative (numeric & variable) metabolite columns
met_names <- met_names[
  vapply(Data1[met_names], function(x) is.numeric(x) && sd(x, na.rm = TRUE) > 0, logical(1))
]

message(sprintf("Will test %d metabolites in adjusted models.", length(met_names)))

# ----------------------------
# Adjusted conditional logistic regression (per metabolite)
# ----------------------------
fit_one_metabolite <- function(met_name) {
  df <- Data1 %>%
    mutate(met = .data[[met_name]]) %>%
    filter(is.finite(met)) %>%
    group_by(indikatorvariabel) %>%
    filter(n() == 2) %>%                # keep complete pairs
    ungroup()
  
  # If nothing left, return NA row
  if (nrow(df) == 0) {
    return(tibble(
      metID   = met_name,
      OR      = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      p       = NA_real_, n_pairs = NA_integer_
    ))
  }
  
  fit <- try(
    clogit(CHD ~ scale(met) + KOEN_BARN + BMI + GESTATIONSALDER_DAGE +
             RYG_cat + alcohol + parity + strata(indikatorvariabel),
           data = df),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(tibble(
      metID   = met_name,
      OR      = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      p       = NA_real_, n_pairs = dplyr::n_distinct(df$indikatorvariabel)
    ))
  }
  
  co <- summary(fit)$coefficients
  if (!"scale(met)" %in% rownames(co)) {
    return(tibble(
      metID = met_name,
      OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      p = NA_real_, n_pairs = dplyr::n_distinct(df$indikatorvariabel)
    ))
  }
  
  b  <- co["scale(met)","coef"]
  se <- co["scale(met)","se(coef)"]
  p  <- co["scale(met)","Pr(>|z|)"]
  
  tibble(
    metID   = met_name,
    OR      = exp(b),
    CI_low  = exp(b - 1.96*se),
    CI_high = exp(b + 1.96*se),
    p       = p,
    n_pairs = dplyr::n_distinct(df$indikatorvariabel)
  )
}

adj_results <- purrr::map_dfr(met_names, fit_one_metabolite) %>%
  mutate(FDR = p.adjust(p, method = "BH"))

# Save adjusted case-control results
write_csv(adj_results, file.path(out_dir, "clogit_adjusted_per_metabolite.csv"))

# Label significance, attach names, list nominal hits

adj_results_labeled <- adj_results %>%
  left_join(annotations %>% dplyr::select(metID, curated_annotation), by = "metID") %>%
  mutate(
    sig_cat = case_when(
      is.na(p)               ~ "Not estimated",
      !is.na(FDR) & FDR < .05 ~ "Statistically significant (FDR < 0.05)",
      p < .05                ~ "Nominally significant (p < 0.05)",
      TRUE                   ~ "Not significant"
    ),
    # Prefer curated name, fall back to metID
    display_name = ifelse(!is.na(curated_annotation) & curated_annotation != "",
                          curated_annotation, metID)
  ) %>%
  arrange(p)

# Quick overview counts
adj_results_labeled %>%
  count(sig_cat) %>%
  arrange(desc(n)) %>%
  print(n = Inf)

# List of nominally significant (p < 0.05 but FDR >= 0.05)
nominal_list <- adj_results_labeled %>%
  filter(sig_cat == "Nominally significant (p < 0.05)") %>%
  pull(metID) %>%
  unique()


cat("\nNominally significant metabolites (CLR, p < 0.05):\n",
    paste(nominal_list, collapse = ", "), "\n\n")


# write labeled table and nominal list
readr::write_csv(adj_results_labeled, file.path(out_dir, "clogit_adjusted_labeled.csv"))
readr::write_lines(nominal_list, file.path(out_dir, "clogit_nominal_metabolites.txt"))


# ----------------------------
# Interaction analysis: metabolite × CHD_subtype
# Notes:
# - Model: CHD ~ scale(met)*CHD_subtype + covariates + strata(pair)
# - LRT comparing interaction vs no interaction (main effects only)
# - Subtype-specific OR per 1-SD for the metabolite (via delta method)
# - Reference subtype = first level of CHD_subtype

# build Data_int with APVR -> Other applied to the pair subtype
Data_int <- Data1 %>%
  group_by(indikatorvariabel) %>%
  mutate(
    .case_subtype = ifelse(CHD == 1, as.character(CHD_type1), NA_character_),
    .pair_subtype = suppressWarnings(dplyr::first(na.omit(.case_subtype)))
  ) %>%
  ungroup() %>%
  filter(!is.na(.pair_subtype)) %>%
  mutate(
    # collapse APVR into Other
    CHD_subtype = case_when(.pair_subtype == "APVR" ~ "Other",
                            .pair_subtype == "Oth_spec_de" ~ "Other",
                            TRUE ~ .pair_subtype),
    CHD_subtype = factor(CHD_subtype)
  ) %>%
  select(-.case_subtype, -.pair_subtype) %>%
  droplevels()


# pick reference level - picked the one with highest N
if ("Septal_def" %in% levels(Data_int$CHD_subtype)) {
  Data_int$CHD_subtype <- relevel(Data_int$CHD_subtype, ref = "Septal_def")
}

# Sanity checks
table(Data_int$CHD, useNA="ifany")
table(Data_int$CHD_subtype, useNA="ifany")
length(levels(Data_int$CHD_subtype))  # must be >= 2


# helper to enforce identical rows for m0/m1 (prevents LRT=NA)
complete_for_interaction <- function(df) {
  # Variables that appear in the *interaction* model
  vars_needed <- c("CHD","KOEN_BARN","BMI","GESTATIONSALDER_DAGE",
                   "RYG_cat","alcohol","parity","indikatorvariabel","CHD_subtype","met")
  df %>%
    mutate(met = as.numeric(met), met_z = as.numeric(scale(met))) %>%
    filter(is.finite(met_z)) %>%
    # complete cases for *all* vars that appear in m1
    filter(stats::complete.cases(dplyr::across(all_of(vars_needed))))
}

fit_interaction_one <- function(met_name) {
  # Start from matched pairs only
  df0 <- Data_int %>%
    mutate(met = .data[[met_name]]) %>%
    group_by(indikatorvariabel) %>%
    filter(n() == 2) %>%
    ungroup()
  
  # Same complete-case rows for both models
  df <- complete_for_interaction(df0)
  
  fac_vars <- c("CHD_subtype","KOEN_BARN","RYG_cat","alcohol","parity")
  fac_ok <- fac_vars[sapply(
    fac_vars,
    function(v) v %in% names(df) && is.factor(df[[v]]) && nlevels(df[[v]]) >= 2
  )]
  
  cont_vars <- c("BMI","GESTATIONSALDER_DAGE")
  cont_ok <- cont_vars[cont_vars %in% names(df)]
  
  # We *must* have CHD_subtype for interaction
  if (!("CHD_subtype" %in% fac_ok)) {
    return(list(
      header = tibble(
        metID = met_name,
        p_met_main = NA_real_,
        p_LRT_interaction = NA_real_,
        n_pairs = dplyr::n_distinct(df$indikatorvariabel)
      ),
      by_subtype = tibble(
        metID = met_name,
        subtype = NA_character_,
        OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        p_subtype = NA_real_
      )
    ))
  }
  
  rhs_m0 <- paste(c("met_z", fac_ok, cont_ok), collapse = " + ")
  rhs_m1 <- paste(c("met_z * CHD_subtype",
                    setdiff(fac_ok, "CHD_subtype"),
                    cont_ok), collapse = " + ")
  
  form_m0 <- as.formula(paste("CHD ~", rhs_m0, "+ strata(indikatorvariabel)"))
  form_m1 <- as.formula(paste("CHD ~", rhs_m1, "+ strata(indikatorvariabel)"))
  
  m0 <- try(clogit(form_m0, data = df), silent = TRUE)
  m1 <- try(clogit(form_m1, data = df), silent = TRUE)
  
  
  if (inherits(m0, "try-error") || inherits(m1, "try-error")) {
    return(list(
      header = tibble(
        metID = met_name,
        p_met_main = NA_real_,
        p_LRT_interaction = NA_real_,
        n_pairs = dplyr::n_distinct(df$indikatorvariabel)
      ),
      by_subtype = tibble(
        metID = met_name,
        subtype = NA_character_,
        OR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        p_subtype = NA_real_
      )
    ))
  }
  
  ## --- LRT from log-likelihoods ---
  ll0 <- as.numeric(logLik(m0))
  ll1 <- as.numeric(logLik(m1))
  df_diff <- attr(logLik(m1), "df") - attr(logLik(m0), "df")
  p_LRT <- if (is.finite(ll0) && is.finite(ll1) && df_diff > 0) {
    pchisq(2 * (ll1 - ll0), df = df_diff, lower.tail = FALSE)
  } else {
    NA_real_
  }
  
  co <- summary(m1)$coefficients
  V  <- vcov(m1)
  
  # p-value for main metabolite effect in interaction model
  p_met_main <- if ("met_z" %in% rownames(co)) {
    co["met_z", "Pr(>|z|)"]
  } else {
    NA_real_
  }
  
  ref_level <- levels(df$CHD_subtype)[1]
  
  ## --- subtype-specific effects: beta and p_subtype ---
  
  # Reference subtype: beta_met only
  out_ref <- if ("met_z" %in% rownames(co)) {
    beta <- co["met_z", "coef"]
    se   <- co["met_z", "se(coef)"]
    z    <- beta / se
    p_ref <- 2 * pnorm(abs(z), lower.tail = FALSE)
    tibble(
      subtype  = ref_level,
      OR       = exp(beta),
      CI_low   = exp(beta - 1.96 * se),
      CI_high  = exp(beta + 1.96 * se),
      p_subtype = p_ref
    )
  } else {
    tibble(
      subtype  = ref_level,
      OR       = NA_real_,
      CI_low   = NA_real_,
      CI_high  = NA_real_,
      p_subtype = NA_real_
    )
  }
  
  # Other subtypes: beta_met + beta_interaction
  others <- setdiff(levels(df$CHD_subtype), ref_level)
  
  out_others <- purrr::map_dfr(others, function(s) {
    term <- paste0("met_z:CHD_subtype", s)
    if (!("met_z" %in% rownames(co)) || !(term %in% rownames(co))) {
      return(tibble(
        subtype  = s,
        OR       = NA_real_,
        CI_low   = NA_real_,
        CI_high  = NA_real_,
        p_subtype = NA_real_
      ))
    }
    b  <- co["met_z", "coef"] + co[term, "coef"]
    var_b <- V["met_z", "met_z"] + V[term, term] + 2 * V["met_z", term]
    se <- sqrt(var_b)
    z  <- b / se
    p  <- 2 * pnorm(abs(z), lower.tail = FALSE)
    tibble(
      subtype  = s,
      OR       = exp(b),
      CI_low   = exp(b - 1.96 * se),
      CI_high  = exp(b + 1.96 * se),
      p_subtype = p
    )
  })
  
  list(
    header = tibble(
      metID = met_name,
      p_met_main = p_met_main,
      p_LRT_interaction = p_LRT,
      n_pairs = dplyr::n_distinct(df$indikatorvariabel)
    ),
    by_subtype = bind_rows(out_ref, out_others) %>%
      mutate(metID = met_name, .before = 1)
  )
}


# Check that every kept pair has a non-missing subtype shared by both rows
Data_int %>%
  group_by(indikatorvariabel) %>%
  summarize(n = n(),
            n_cases = sum(CHD == 1),
            n_ctrls = sum(CHD == 0),
            n_subtypes = n_distinct(CHD_subtype),
            subtype_example = dplyr::first(CHD_subtype)) %>%
  count(n, n_cases, n_ctrls, n_subtypes)

# Expect: n = 2, n_cases = 1, n_ctrls = 1, n_subtypes = 1 

# Sanity-check fit interaction on one metabolite end-to-end
one <- met_names[1]
tmp <- fit_interaction_one(one)
tmp$header
tmp$by_subtype

safe_fit <- purrr::safely(fit_interaction_one, otherwise = NULL, quiet = TRUE)
res_list <- purrr::map(met_names, safe_fit)

# Headers
interaction_header <- res_list %>%
  keep(~ !is.null(.x$result)) %>%
  purrr::map("result") %>%
  purrr::map("header") %>%
  bind_rows() %>%
  mutate(
    FDR_LRT      = p.adjust(p_LRT_interaction, method = "BH"),
    FDR_met_main = p.adjust(p_met_main,         method = "BH"),
    sig_label = case_when(
      is.na(p_LRT_interaction)          ~ "No test",
      FDR_LRT < 0.05                    ~ "Significant interaction (FDR < 0.05)",
      p_LRT_interaction < 0.05          ~ "Nominal interaction (p < 0.05)",
      TRUE                              ~ "NS"
    )
  )

# By-subtype ORs + p-values
interaction_by_subtype <- res_list %>%
  keep(~ !is.null(.x$result)) %>%
  purrr::map("result") %>%
  purrr::map("by_subtype") %>%
  bind_rows() %>%
  group_by(subtype) %>%
  mutate(FDR_subtype = p.adjust(p_subtype, method = "BH")) %>%
  ungroup()


message("\n--- Adjusted case-control (clogit) summary ---")
print(adj_results %>%
        transmute(metID, OR_CI = sprintf("%.2f [%.2f, %.2f]", OR, CI_low, CI_high),
                  p = signif(p, 3), FDR = signif(FDR, 3)) %>%
        arrange(p) %>% head(10))

message("\n--- Top interaction LRT (lowest p) ---")
print(interaction_header %>% arrange(p_LRT_interaction) %>% head(100))
