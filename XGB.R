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
gwas_data_rf <- read_delim("Selected_SNP.txt")
#gwas_data_rf <- read_delim("Selected_SNP_PH.txt")
# --- Regression Task --------------------------------------------------------
set.seed(123)
task_gwas_rf <- TaskRegr$new(
  id = "gwas_rf",
  backend = gwas_data_rf,
  target = "TKW" # change to "PH" 
)

# ============================================================
# 1. XGBoost Learner & Hyperparameter Tuning
# ============================================================

# Define learner and hyperparameter search space
learner_rf <- lrn("regr.xgboost",
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
tuner_rf <- tnr("mbo")
measure_mse <- msr("regr.mse")
terminator_eval <- trm("evals", n_evals = 10)

# Define resampling strategies
outer_rs <- rsmp("cv", folds = 5)  # outer loop
inner_rs <- rsmp("cv", folds = 5)  # inner loop

# Nested CV tuner
auto_tuner_rf <- auto_tuner(
  tuner = tuner_rf,
  learner = learner_rf,
  resampling = inner_rs,
  measure = measure_mse,
  terminator = terminator_eval
)

# --- Run Nested Cross-Validation --------------------------------------------
rr <- resample(task_gwas_rf, auto_tuner_rf, outer_rs, store_models = TRUE)

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
learner_rf$param_set$values <- c(best_hyperpars[[1]], list(nrounds = 1000))

# ============================================================
# 3. Variable Importance (Average over Multiple Runs)
# ============================================================

n_runs <- 100
feature_names <- task_gwas_rf$feature_names
importance_list <- vector("list", n_runs)

for (i in seq_len(n_runs)) {
  learner_rf$train(task_gwas_rf)
  imp_matrix <- xgb.importance(feature_names = feature_names, model = learner_rf$model)
  importance_list[[i]] <- as.data.frame(imp_matrix)
  cat("Run", i, "completed\n")
}

# Average importance scores
avg_importance <- bind_rows(importance_list) %>%
  group_by(Feature) %>%
  summarize(
    Gain = mean(Gain, na.rm = TRUE),
    Cover = mean(Cover, na.rm = TRUE),
    Frequency = mean(Frequency, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  rename(SNP = Feature)

# ============================================================
# 4. Permutation Test for Empirical Threshold
# ============================================================

n_perm <- 1000
perm_max_importances <- numeric(n_perm)
feature_names <- setdiff(task_gwas_rf$feature_names, task_gwas_rf$target_names)
learner_rf_perm <- learner_rf$clone(deep = TRUE)

for (i in seq_len(n_perm)) {
  data_perm <- as.data.table(task_gwas_rf$data())
  data_perm[, TKW := sample(TKW)]  # change to PH if target is PH
  
  task_perm <- TaskRegr$new(
    id = paste0("gwas_rf_perm_", i),
    backend = data_perm,
    target = "TKW" # change to PH if needed
  )
  
  learner_rf_perm$train(task_perm)
  imp_perm <- xgb.importance(feature_names = feature_names, model = learner_rf_perm$model)
  perm_max_importances[i] <- max(imp_perm$Gain)
  cat("Permutation", i, "completed\n")
}

# 95th percentile empirical threshold
emp_threshold_perm <- quantile(perm_max_importances, probs = 0.95)
cat("Empirical 95th percentile threshold:", emp_threshold_perm, "\n")

# ============================================================
# 5. Merge with Marker Metadata and Save Results
# ============================================================

#yGM <- read.csv("myGM.csv")

#importance_df2 <- merge(avg_importance, myGM, by = "SNP")
#write.table(importance_df2, "importance_xgboost_PH.txt", row.names = FALSE)
#write.table(emp_threshold_perm, "emp_threshold_xgboost_TKW.txt", row.names = FALSE)

# Optional: save workspace
# save.image("mlr3_xgboost_PH.RData")
