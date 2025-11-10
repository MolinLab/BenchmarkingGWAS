# ============================================================
# Ranger Random Forest GWAS Model with Nested CV, Variable Importance & Permutation Test
# ============================================================

# --- Setup ------------------------------------------------------------------
setwd()

# --- Load Required Packages -------------------------------------------------
library(mlr3verse)
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(paradox)
library(mlr3misc)
library(mlr3mbo)
library(mlr3viz)
library(data.table)
library(readr)
library(dplyr)
library(ranger)
library(tidyr)

# --- Data -------------------------------------------------------------------
# Selected SNPs by GLMnet
gwas_data_rf <- read_delim("Selected_SNP.txt")
# gwas_data_rf <- read_delim("Selected_SNP_PH.txt")

# --- Regression Task --------------------------------------------------------
set.seed(123)
task_gwas_rf <- TaskRegr$new(
  id = "gwas_rf",
  backend = gwas_data_rf,
  target = "TKW" # change to "PH" or trait of interest
)

# ============================================================
# 1. Ranger Learner & Hyperparameter Tuning
# ============================================================

learner_rf <- lrn(
  "regr.ranger",
  num.threads = 1,
  importance = "impurity",
  num.trees = to_tune(300, 1500),
  mtry = to_tune(1, 300),
  min.node.size = to_tune(1, 500),
  max.depth = to_tune(1, 100)
)

# Use Bayesian optimization
tuner_rf <- tnr("mbo")
measure_rmse <- msr("regr.rmse")

# Define resampling strategies
outer_rs <- rsmp("cv", folds = 5)  # outer loop
inner_rs <- rsmp("cv", folds = 5)  # inner loop

# Tuning control
terminator_eval <- trm("evals", n_evals = 10)

# Nested CV auto-tuner
auto_tuner_rf <- auto_tuner(
  tuner = tuner_rf,
  learner = learner_rf,
  resampling = inner_rs,
  measure = measure_rmse,
  terminator = terminator_eval
)

# --- Run Nested Cross-Validation --------------------------------------------
rr <- resample(task_gwas_rf, auto_tuner_rf, outer_rs, store_models = TRUE)
autoplot(rr)

# ============================================================
# 2. Model Performance Evaluation
# ============================================================

measures <- msrs(c(
  "regr.rmse", "regr.rsq", "regr.mae", "regr.pbias",
  "time_train", "time_predict", "time_both"
))

# Aggregate results
rr$aggregate(measures)

# Identify best-performing outer fold
score_dt <- rr$score(msr("regr.mse"))
best_fold_index <- which.min(score_dt$regr.mse)

cat("Best outer fold index:", best_fold_index, "\n")

# Extract best hyperparameters
best_hyperpars <- rr$learners[[best_fold_index]]$tuning_result$x
cat("Best hyperparameters from the best outer fold:\n")
print(best_hyperpars)

# Update learner with best hyperparameters
learner_rf$param_set$values <- c(best_hyperpars[[1]], list(importance = "impurity"))

# ============================================================
# 3. Variable Importance (Average over Multiple Runs)
# ============================================================

n_runs <- 100
feature_names <- task_gwas_rf$feature_names
importance_list <- vector("list", n_runs)

for (i in seq_len(n_runs)) {
  learner_rf$train(task_gwas_rf)
  
  importance_scores <- learner_rf$model$variable.importance
  importance_df <- data.frame(
    Feature = names(importance_scores),
    Importance = as.vector(importance_scores)
  )
  
  importance_list[[i]] <- importance_df
  cat("Run", i, "completed\n")
}

# Average importance scores
avg_importance <- bind_rows(importance_list) %>%
  group_by(Feature) %>%
  summarize(Importance = mean(Importance, na.rm = TRUE)) %>%
  ungroup() %>%
  rename(SNP = Feature)
# ============================================================
# 4. Permutation Test for Empirical Threshold
# ============================================================

n_perm <- 1000
perm_max_importances <- numeric(n_perm)
learner_rf_perm <- learner_rf$clone(deep = TRUE)

for (i in seq_len(n_perm)) {
  data_perm <- as.data.table(task_gwas_rf$data())
  
  # Permute the target variable
  data_perm[, TKW := sample(TKW)]  # change to PH 
  
  # Create new permuted task
  task_perm <- TaskRegr$new(
    id = paste0("gwas_rf_perm_", i),
    backend = data_perm,
    target = "TKW" # change to PH 
  )
  
  learner_rf_perm$train(task_perm)
  
  perm_max_importances[i] <- max(learner_rf_perm$model$variable.importance, na.rm = TRUE)
  cat("Permutation", i, "completed\n")
}

# 95th percentile empirical threshold
emp_threshold_perm <- quantile(perm_max_importances, probs = 0.95)
cat("Empirical 95th percentile threshold (max importances):", emp_threshold_perm, "\n")

# Identify significant features
significant_features <- avg_importance$SNP[avg_importance$Importance > emp_threshold_perm]
cat("Number of features exceeding the empirical threshold:", length(significant_features), "\n")



# ============================================================
# 5. Save Results
# ============================================================

#save.image("mlr3_RF_TKW.RData")

# Optionally merge with marker metadata
# myGM <- read.csv("../../../imputation_test/new_data/myGM_PH.csv")
# importance_df2 <- merge(imp_summary, myGM, by = "SNP")
# write.table(importance_df2, "importance_RF_PH.txt", row.names = FALSE)
# write.table(imp_threshold, "emp_threshold_RF_PH.txt", row.names = FALSE)

 # ============================================================
# average importance across folds
# ============================================================
n_reps <- 100                # how many repeated 5-fold CV runs
n_workers_reps <- 10             # how many parallel workers for the 100 reps 

feature_names <- task_gwas_rf$feature_names
n_features <- length(feature_names)

run_one_rep <- function(rep_idx) {
  cat("Starting repetition", rep_idx, "\n")
  set.seed(1000 + rep_idx)
  
  learner_clone <- learner_rf$clone(deep = TRUE)
  resampling_clone <- outer_rs$clone(deep = TRUE)
  
  rr_local <- resample(
    task = task_gwas_rf,
    learner = learner_clone,
    resampling = resampling_clone,
    store_models = TRUE
  )
  
  per_fold_importances <- matrix(
    0, nrow = n_features, ncol = length(rr_local$learners),
    dimnames = list(feature_names, paste0("fold", seq_along(rr_local$learners)))
  )
  
  for (i in seq_along(rr_local$learners)) {
    learner_i <- rr_local$learners[[i]]
    imp_vec <- learner_i$importance()
    
    tmp <- numeric(n_features)
    names(tmp) <- feature_names
    
    if (!is.null(imp_vec) && length(imp_vec) > 0)
      tmp[names(imp_vec)] <- as.numeric(imp_vec)
    
    per_fold_importances[, i] <- tmp
  }
  
  mean_imp_rep <- rowMeans(per_fold_importances, na.rm = TRUE)
  cat("Finished repetition", rep_idx, "\n")
  
  return(list(mean_importance = mean_imp_rep))
}

plan(multisession, workers = n_workers_reps)  # parallel across reps
rep_results <- future_lapply(seq_len(n_reps), function(r) {
  run_one_rep(r)
}, future.seed = TRUE)
plan(sequential)

imp_all_runs <- sapply(rep_results, function(res) res$mean_importance)


imp_mean <- rowMeans(imp_all_runs, na.rm = TRUE)

imp_summary <- data.table(
  marker = feature_names,
  mean_importance = imp_mean)
                       [order(-mean_importance)]




B <- 1000 #number permutations
n_workers <- 10

compute_mean_imp_from_rr <- function(rr_local, feature_names) {
  n_iter_local <- length(rr_local$learners)
  imp_mat_local <- matrix(0, nrow = length(feature_names), ncol = n_iter_local,
                          dimnames = list(feature_names, paste0("iter", seq_len(n_iter_local))))
  for (i in seq_len(n_iter_local)) {
    imp_vec <- rr_local$learners[[i]]$importance()
    if (!is.null(imp_vec) && length(imp_vec) > 0) {
      tmp <- numeric(length(feature_names))
      names(tmp) <- feature_names
      tmp[names(imp_vec)] <- as.numeric(imp_vec)
      imp_mat_local[, i] <- tmp
    }
  }
  rowMeans(imp_mat_local, na.rm = TRUE)
}

run_one_perm <- function(seed = NULL, resampling_clone) {
  if (!is.null(seed)) set.seed(seed)
  
  data_perm <- as.data.table(task_gwas_rf$data())
  target_name <- task_gwas_rf$target_names
  data_perm[[target_name]] <- sample(data_perm[[target_name]])
  
  task_perm <- TaskRegr$new(
    id = paste0("perm_task_", sample.int(1e8, 1)),
    backend = data_perm,
    target = target_name
  )
  
  learner_clone <- learner_rf$clone(deep = TRUE)
  rr_perm <- resample(
    task = task_perm,
    learner = learner_clone,
    resampling = resampling_clone,
    store_models = TRUE
  )
  
  mean_imp_perm <- compute_mean_imp_from_rr(rr_perm, feature_names)
  max(mean_imp_perm, na.rm = TRUE)
}

set.seed(123)
seeds <- sample.int(.Machine$integer.max, B)
resampling_clone <- outer_rs$clone(deep = TRUE)

cat("Starting permutation threshold with", B, "permutations on", n_workers, "workers\n")

future::plan("multisession", workers = n_workers)
with_progress({
  perm_max_imp <- future_lapply(seq_len(B), function(b) {
    run_one_perm(seed = seeds[b], resampling_clone = resampling_clone)
  }, future.seed = TRUE)
})
future::plan("sequential")

perm_max_imp <- vapply(perm_max_imp, as.numeric, numeric(1))
imp_threshold <- as.numeric(quantile(perm_max_imp, probs = alpha, na.rm = TRUE))

cat(sprintf("permutation threshold (%.1f%% percentile): %g\n",0.95, imp_threshold))
