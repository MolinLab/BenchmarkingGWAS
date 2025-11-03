# ============================================================
# XGBoost GWAS Model with Nested CV, Variable Importance & Permutation Test
# ============================================================

# --- Setup ------------------------------------------------------------------
setwd()

# --- Load Required Packages -------------------------------------------------
library(mlr3verse)
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3misc)
library(mlr3mbo)
library(mlr3viz)
library(paradox)
library(data.table)
library(xgboost)
library(readr)
library(dplyr)

# --- Data -------------------------------------------------------------------
# Selected SNPs by GLMnet
gwas_data_xgb <- read_delim("Selected_SNP.txt")
#gwas_data_xgb <- read_delim("Selected_SNP_PH.txt")
# --- Regression Task --------------------------------------------------------
set.seed(123)
task_gwas_xgb <- TaskRegr$new(
  id = "gwas_xgb",
  backend = gwas_data_xgb,
  target = "TKW" # change to "PH" 
)

# ============================================================
# 1. XGBoost Learner & Hyperparameter Tuning
# ============================================================

# Define learner and hyperparameter search space
learner_xgb <- lrn("regr.xgboost",
                  nthread = 1,
                  eta = to_tune(0.0001, 0.3),# slow done learning rate to prevent overfitting
                  gamma = to_tune(0, 5),# minimum amount of splitting by which a node must improve the loss function
                  max_depth = to_tune(1, 10),# maximum number of levels deep that each tree can grow
                  colsample_bytree = to_tune(0.1, 1),#proportion of predictor variables sampled for each tree
                  min_child_weight = to_tune(0, 10),#minimum degree of impurity needed in a node before attempting to split it
                  lambda = to_tune(1e-3, 1000),#L2 regularization term on weights. Increasing this value will make model more conservative.
                  alpha = to_tune(1e-3, 1000),#L1 regularization term on weights. Increasing this value will make model more conservative.
                  nrounds = 1000# number of sequentially built trees
)

# Define tuner, measure, and terminator
tuner_xgb <- tnr("mbo")
measure_mse <- msr("regr.mse")
terminator_eval <- trm("evals", n_evals = 10)

# Define resampling strategies
outer_rs <- rsmp("cv", folds = 5)  # outer loop
inner_rs <- rsmp("cv", folds = 5)  # inner loop

# Nested CV tuner
auto_tuner_xgb <- auto_tuner(
  tuner = tuner_xgb,
  learner = learner_xgb,
  resampling = inner_rs,
  measure = measure_mse,
  terminator = terminator_eval
)

# --- Run Nested Cross-Validation --------------------------------------------
rr <- resample(task_gwas_xgb, auto_tuner_xgb, outer_rs, store_models = TRUE)

# ============================================================
# 2. Model Performance Evaluation
# ============================================================

# Define performance measures
measures <- msrs(c(
  "regr.rmse", "regr.mae", "regr.rsq",
  "regr.pbias", "time_train", "time_predict", "time_both"
))

# Aggregate results
rr(aggregate(measures = measures))

# Identify best-performing outer fold
score_dt <- rr$score(msr("regr.mse"))
best_fold_index <- which.min(score_dt$regr.mse)

cat("Best outer fold index:", best_fold_index, "\n")

# Extract best hyperparameters
best_hyperpars <- rr$learners[[best_fold_index]]$tuning_result$x
cat("Best hyperparameters from the best outer fold:\n")
print(best_hyperpars)

# Update learner with best hyperparameters
learner_xgb$param_set$values <- c(best_hyperpars[[1]], list(nrounds = 1000))

# ============================================================
# 3. Variable Importance (Average over Multiple Runs)
# ============================================================

n_reps <- 100                # how many repeated 5-fold CV runs
n_workers_reps <- 10              # how many parallel workers for the 100 reps (adjust to your resources)

feature_names <- GY_genotype_task$feature_names
n_features <- length(feature_names)

run_one_rep <- function(rep_idx) {
  cat("Starting repetition", rep_idx, "\n")
  set.seed(1000 + rep_idx)
  
  learner_clone <- learner_xgb$clone(deep = TRUE)
  resampling_clone <- outer_rs$clone(deep = TRUE)
  
  rr_local <- resample(
    task = tasg_gwas_xgb,
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
    imp_vec <- learner_i$importance() # here we used xgb.importance function. Using learner_i$importance we get an easier way of extracting gain 
    
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
  SNP = feature_names,
  mean_importance = imp_mean
  )[order(-mean_importance)]


# ============================================================
# 4. Permutation Test for Empirical Threshold
# ============================================================
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
  
  data_perm <- as.data.table(task_gwas_xgb$data())
  target_name <- task_gwas_xgb$target_names
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
emp_threshold_perm <- as.numeric(quantile(perm_max_imp, probs = alpha, na.rm = TRUE))

cat(sprintf("permutation threshold (%.1f%% percentile): %g\n",0.95, emp_threshold_perm))
# ============================================================
# 5. Merge with Marker Metadata and Save Results
# ============================================================

#yGM <- read.csv("myGM.csv")

#importance_df2 <- merge(imp_summary, myGM, by = "SNP")
#write.table(importance_df2, "importance_xgboost_PH.txt", row.names = FALSE)
#write.table(emp_threshold_perm, "emp_threshold_xgboost_TKW.txt", row.names = FALSE)

# Optional: save workspace
# save.image("mlr3_xgboost_PH.RData")
