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
# importance_df2 <- merge(avg_importance, myGM, by = "SNP")
# write.table(importance_df2, "importance_RF_PH.txt", row.names = FALSE)
# write.table(emp_threshold_perm, "emp_threshold_RF_PH.txt", row.names = FALSE)
