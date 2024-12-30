# Author: Luca Pestarino
# Created: 24.04.2023
# Updated: 30.12.2024
# Main script used to validate Umu et al., 2022 paper results using only selected biomarkers

rm(list=ls())

# Install needed packages
packages <- c("xgboost", "future", "pROC", "caret", "plotROC", 
              "optmatch", "dplyr", "purrr", "tibble", "stringr", 
              "tidyr", "tidyverse", "mlr", "readxl")
installed <- packages %in% installed.packages()
if (any(!installed)) {
  install.packages(packages[!installed])
}


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!"DESeq2" %in% installed.packages()) {
  BiocManager::install("DESeq2")
}

lapply(c(packages, "DESeq2"), library, character.only = TRUE)

# Load functions file
source("~/LC_biomarkers_validation_functions.R") # Path to the functions file

# Data files
data_dir <- "~/Data/" # directory with data files
count_tables_dir <- "~/count_tables/" # directory with the count tables only
save_dir <- "~/Results/" # directory where to save the results


caco_file <- readRDS("~/lc_main_diff_condition_with_nonsmokers_all_lc_df.rds") # File with clinical information
smokers_caco_file <- readRDS("~/lc_main_clinical_smoking_updated_df.rds") # File with clinical information for smokers only (used for the training process)
non_smokers_caco_file <- caco_file %>% filter(!sample %in% smokers_caco_file$sample) # Extract non-smokers samples

# Print median age at diagnosis
print(median((caco_file$TDATE_DIAG_TIME))/52.1429)


#####
# Apply women filter
#caco_file <- caco_file %>% filter(sex == "F") %>% filter(condition == "LC")
#smokers_caco_file <- smokers_caco_file %>% filter(sex == "F")
#matched_variables <- c("BDg", "age_continuous") #### list of variables to match on
# non_smokers_caco_file <- caco_file %>% filter(!sample %in% smokers_caco_file$sample) # Extract non-smokers for second testing

#print(median((caco_file$TDATE_DIAG_TIME))/52.1429)



# Load biomarkers lists from original LC paper
mirna_bm <- read_excel(paste(data_dir, "Supplemental_Table_1.xlsx", sep = ""), sheet = "miRNA models")
miscrna_bm <- read_excel(paste(data_dir, "Supplemental_Table_1.xlsx", sep = ""), sheet = "miscRNA models")
all_bm <- read_excel(paste(data_dir, "Supplemental_Table_1.xlsx", sep = ""), sheet = "All RNAs models")

bm_list <- list(miRNA = mirna_bm, misc_RNA = miscrna_bm, All_RNA = all_bm) 

# Important variables to define
matched_variables <- c("BDg", "age_continuous", "sex") #### list of variables to match the samples on
target_variable <- "condition"   #### name of the diagnosis condition variable: LC vs C
smokers_caco_file <- smokers_caco_file %>% filter((diagnosetimeinyears < 10 & diagnosetimeinyears >= 0) | timetodiagnose %in% "C") %>%
  mutate_at(c("sex", "BDg"), as.factor) %>% as.data.frame() #### select only pre-diagnostic samples 0-10 years before diagnosis
cancer_type <- c("Control", "NSCLC", "Others", "SCLC") #### check if the naming for different cancer types is compatible

seeds <- 1:5


#####
# load count tables or create them if not existing 
####

if (file.exists(paste(data_dir, "rna_list_end_to_end.rds", sep = ""))) {
  rna_df_list <- readRDS(file = paste(data_dir, "rna_list_end_to_end.rds", sep = ""))
} else {
  files_list <- list.files(count_tables_dir)
  rna_df_list <- list()
  rna_list <- lapply(strsplit(files_list, split="[.]"), '[[', 1)
  rna_list <- gsub("protein_coding", "mRNA_fragments", rna_list)
  for (rna_file in files_list){
    file_name <- strsplit(rna_file, split="[.]")[[1]][1]
    count_table <- read.delim(file = paste(count_tables_dir, rna_file, sep = ""))
    rna_df_list[[file_name]] <- count_table
  }
  names(rna_df_list) <- rna_list
  saveRDS(rna_df_list, paste(data_dir, "rna_list_end_to_end.rds", sep = ""))
}


# Normalise the count tables
if (file.exists(paste(data_dir, "norm_rna_list_no_filter_end_to_end.rds", sep = ""))){
  norm_rna_list <- readRDS(file = paste(data_dir, "norm_rna_list_no_filter_end_to_end.rds", sep = ""))
} else {
  norm_rna_list <- lapply(rna_df_list, deseq2_normalise, caco_file, threshold = 5, filtering = F)
  saveRDS(norm_rna_list, paste(data_dir, "norm_rna_list_no_filter_end_to_end.rds", sep = ""))
}

norm_rna_list$gencode <- NULL # Remove gencode class to avoid duplicate RNA IDs


#####
### Main Analyses

parameters_grid <- makeParamSet(makeDiscreteParam("subsample", values = c(0.3, 0.4, 0.55, 0.75)),
                                makeDiscreteParam("min_child_weight", values = c(1, 2, 6)),
                                makeDiscreteParam("colsample_bytree", values = c(0.6, 1)),
                                makeDiscreteParam("eta", values = c(0.05, 0.1, 0.3)),
                                makeDiscreteParam("max_depth", values = c(3, 6, 10)),
                                makeDiscreteParam("early_stopping_rounds", values = c(700)),
                                makeDiscreteParam("gamma", values = c(1, 5, 10)),
                                makeDiscreteParam("lambda", values = c(1)),
                                makeDiscreteParam("alpha", values = c(1, 10)),
                                makeDiscreteParam("nrounds", values = c(1000)))


final_results_smokers <- lapply(c("All histologies", "NSCLC", "SCLC"), get_models_results, rna_classes = names(bm_list), rna_bm_list = bm_list, meta_file = smokers_caco_file, 
                                target_variable = target_variable, normalised_rna_list = norm_rna_list, matched_variables = matched_variables, 
                                parameters_grid = parameters_grid, save_directory = save_dir, matched_data = F) %>% 
  setNames(c("All histologies", "NSCLC", "SCLC"))

final_results_smokers_by_rna <- purrr::transpose(final_results_smokers)
save(final_results_smokers_by_rna, file = paste(save_dir, "final_results_smokers.RData", sep = ""))


print("Script completed successfully")
