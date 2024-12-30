# Author: Luca Pestarino
# Created: 24.04.2023
# Updated: 30.12.2024
# Functions needed for the script "LC_biomarkers_validation_main.R" to validate Umu et al., 2022 paper results

#####
# Normalise data
deseq2_normalise = function(raw_counts_table, clinical_df, threshold = 5, filtering = T){
  if(filtering == F){threshold = 0}
  colnames(raw_counts_table)[1] <- "ID"
  raw_counts_table$ID <- as.character(raw_counts_table$ID)
  rownames(clinical_df) <- clinical_df$sample #assign row names
  clinical_ids <- clinical_df$sample
  filtered_table <- raw_counts_table %>% group_by(ID) %>% dplyr::select(one_of(clinical_ids)) %>% gather(sample, read.counts, one_of(clinical_ids)) %>% 
    mutate(quant = quantile(read.counts, 0.20)) %>% filter(quant >= threshold) %>% dplyr::select(ID, sample, read.counts) %>% spread(sample, read.counts) 
  filtered_counts_df <- data.frame(filtered_table %>% ungroup() %>% dplyr::select(-ID), row.names = filtered_table$ID)
  if(any(colSums(filtered_counts_df) == 0) | any(rowSums(filtered_counts_df) == 0)){
    filtered_counts_df <- filtered_counts_df + 1
  }
  filtered_counts_df <- filtered_counts_df[, clinical_ids]
  dds <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(countData = filtered_counts_df, colData = clinical_df, design = ~1))
  vst <- DESeq2::varianceStabilizingTransformation(dds, blind = T, fitType = "local")
  return(vst_df = as.data.frame(SummarizedExperiment::assay(vst)))
}

#####
# Plot ROC Curves
####

get_roc_plot <- function(results_df, s = 4, save_dir, plot_name = "", histology_label = "", values_roc = "roc_values", auc_test = "test_auc") {
  df <- results_df %>% select(seed, results) %>% rowwise() %>% mutate(rocfit = list(roc = results[[values_roc]])) %>% select(-results) %>% group_by(seed) %>%
    mutate(n = row_number()) %>% ungroup() %>% hoist(rocfit, Specificity = "specificity", Sensitivity = "sensitivity") %>% unnest(c(Specificity, Sensitivity)) %>%
    select(-rocfit) 
  
  #annotation_df <- results_df %>% summarise(l = t.test(results_df[[auc_test]])$conf.int[1], r = t.test(results_df[[auc_test]])$conf.int[2], m = mean(results_df[[auc_test]]))  %>% 
  #  mutate(ci = paste0("AUC: ", format(m, digits=2)," (95% CI, ", format(l, digits=2), "-", format(r, digits=2), ")"))
  
  unique_values <- unique(results_df[[auc_test]])
  if(length(unique_values) == 1) {
    l <- unique_values
    r <- unique_values
  } else {
    t_test_result <- t.test(results_df[[auc_test]])
    l <- t_test_result$conf.int[1]
    r <- t_test_result$conf.int[2]
  }
  
  annotation_df <- results_df %>%
    summarise(l = l, r = r, m = mean(results_df[[auc_test]])) %>%
    mutate(ci = paste0("AUC: ", format(m, digits = 2), " (95% CI, ", format(l, digits = 2), "-", format(r, digits = 2), ")"))
  
  roc_plot <- ggplot(df, aes(x = 1 - Specificity, y = Sensitivity, color = as.factor(seed))) + 
    geom_line(size = 0.8) + 
    geom_segment(x = 0, xend= 1, y = 0, yend = 1, color = "black", linetype = "dashed") + 
    theme(legend.position = "none", panel.background = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5, size = 14)) +
    ggtitle(histology_label) + 
    geom_text(data = annotation_df, aes(x = 1, y = 0, label = ci, hjust = "inward", vjust = -4), size = s, inherit.aes = FALSE) +
    scale_color_brewer(palette = "Set1")
  
  # pdf(paste(save_dir, plot_name, "_roc_plot.pdf", sep = ""), width = 7, height = 7)
  # print(roc_plot)
  # dev.off()
  
  return(roc_plot)
}


#####
# Train boosting model
####
train_boosting_cv <- function(data, target_variable, seed, tuning_parameters){
  set.seed(seed)
  trainIndex <- createDataPartition(data[, target_variable], p=0.7, times = 1, list=FALSE)
  training_set <-  data[trainIndex,]
  test_set <- data[-trainIndex,]
  
  lrn <- makeLearner("classif.xgboost", predict.type = "prob", predict.threshold = 0.5, par.vals = list(eval_metric = "auc",
                                                                                                        objective = "binary:logistic",
                                                                                                        booster = "gbtree",
                                                                                                        maximize = T,
                                                                                                        nthread = 30))
  
  cv_type <- makeResampleDesc(method = "CV", stratify = T, iters = 5)
  traintask <- makeClassifTask(data = training_set, target = target_variable)
  testtask <- makeClassifTask (data = test_set, target = target_variable)
  
  tuned_model <- tuneParams(learner = lrn, task = traintask, par.set = tuning_parameters, show.info = T, control = makeTuneControlGrid(),
                            resampling = cv_type)
  
  
  lrn_tune <- setHyperPars(lrn, par.vals = tuned_model$x)
  best_model <- train(learner = lrn_tune, task = traintask)
  importances <- as.list(getFeatureImportance(best_model)$res %>% arrange(desc(importance)))
  train_pred <- predict(best_model, traintask)$data$prob.1
  train_y_true <- training_set[, target_variable]
  test_pred <- predict(best_model, testtask)
  performance <- mlr::performance(test_pred, measures = list(auc, acc))
  confusion <- confusionMatrix(test_pred$data$response, test_set[, target_variable])
  probs_pred <- test_pred$data$prob.1
  y_true <- test_set[, target_variable]
  roc_data <- as.list(coords(pROC::roc(response = test_set[, target_variable], predictor = test_pred$data$prob.1, direction = "auto")))
  
  performance_df = tibble(seed = seed, test_auc = performance[[1]], test_acc = performance[[2]]) %>% rowwise() %>% distinct() %>% 
    rowwise() %>% mutate(results = list(list(features_importance = importances, roc_values = roc_data,
                                             probs = probs_pred, test_vals = y_true, best_model = best_model$learner$par.vals)))
  
  return(performance_df)
}


#####
# Wrapper functions
####
get_results <- function(rna_data_list, rna_name, biomarker_list, histology, matched_meta_data, target_variable, parameters_grid = parameters_grid, seeds = seeds,
                        save_directory){
  
  biomarker_list_rna <- biomarker_list[[rna_name]] %>% filter(model == histology)
  
  if(rna_name == "All_RNA" | rna_name == "misc_RNA") {
    rna_data <- do.call(bind_rows, rna_data_list) %>% rownames_to_column("RNA_ID") %>% filter(RNA_ID %in% biomarker_list_rna$Feature) %>% column_to_rownames("RNA_ID")
    
  } else  {
    rna_data <- rna_data_list[[rna_name]] %>% rownames_to_column("RNA_ID") %>% filter(RNA_ID %in% biomarker_list_rna$Feature) %>% column_to_rownames("RNA_ID")
  }
  
  if (all(biomarker_list_rna$Feature %in% rownames(rna_data))) {
    print(paste("All values in", rna_name, histology, "are present in the data."))
  } else {
    missing_rna <- biomarker_list_rna$Feature[!(biomarker_list_rna$Feature %in% rownames(rna_data))]
    stop(paste("Error: missing biomarkers in ", rna_name, histology, "and the missing RNA are: \n", paste(missing_rna, collapse = ", ")))
  }
  
  matched_rna <- rna_data %>% t() %>% data.frame() %>% rownames_to_column("ID") %>% filter(ID %in% matched_meta_data$sample) %>% 
    left_join(matched_meta_data %>% mutate(ID = sample) %>% select(c(ID, condition))) %>% column_to_rownames("ID") %>% mutate(across(c(condition), factor))
  
  performance_results <- lapply(seeds, train_boosting_cv, data = matched_rna, target_variable = "condition", tuning_parameters = parameters_grid) %>% bind_rows()
  
  plot_data <- performance_results %>% select(c(seed, results)) %>% rowwise() %>% mutate(probs= list(roc = results[["probs"]])) %>%
    mutate(test_values = list(roc = results [["test_vals"]])) #%>% unnest(c(probs, test_values))
  
  rna_roc <- get_roc_plot(performance_results, save_dir = save_directory, plot_name = paste(rna_name, "_roc", sep = ""), histology_label = histology)

  return(list(combined_plot = rna_roc, plot_data = plot_data))
}


get_models_results <- function(histology, rna_classes, rna_bm_list, normalised_rna_list, meta_file, target_variable, matched_variables, parameters_grid, 
                               save_directory, matched_data = F){
  if (histology == "NSCLC"){
    cancer_type <- c("Control", "NSCLC")
    save_dir_hist <- paste(save_directory, "NSCLC_", sep = "")
  } else if( histology == "SCLC"){
    cancer_type <- c("Control", "SCLC")
    save_dir_hist <- paste(save_directory, "SCLC_", sep = "")
  } else {
    histology <- "All histologies"
    cancer_type <- c("Control", "NSCLC", "Others", "SCLC")
    save_dir_hist <- paste(save_directory, "All_histologies_", sep = "")
  }
  
  filtered_meta_data <- meta_file %>% filter(cancertype %in% cancer_type)
  
  matched_meta_data <- filtered_meta_data  %>% mutate_at(c(target_variable), as.numeric) %>% mutate(condition = .[, target_variable] -1) #0 is controls, 1 is LC cases
  
  if (matched_data == F) {
    matching_formula <- as.formula(paste0(target_variable, " ~ ", paste(matched_variables, collapse = " + ")))
    matches <- pairmatch(matching_formula, controls = 1, data = matched_meta_data)
    not_matched <- matched_meta_data %>% bind_cols(matches %>% data.frame(groups=.,rows=names(.))) %>% filter(is.na(groups)) %>% 
      dplyr::select(-rows) %>% dplyr::select(c(sample, condition))
    matched_meta_data <- matched_meta_data %>% bind_cols(matches %>% data.frame(groups=.,rows=names(.))) %>% filter(!is.na(groups)) %>% dplyr::select(-rows)
  }
  
  rna_results <- lapply(rna_classes, get_results, rna_data_list = normalised_rna_list, biomarker_list = rna_bm_list, 
                        matched_meta_data = matched_meta_data, target_variable = target_variable, parameters_grid = parameters_grid, seeds = seeds, 
                        save_directory = save_dir_hist, histology = histology) %>% setNames(rna_classes)
  print("done")
  
  return(rna_results)
}
