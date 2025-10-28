# ============================================================
# Elastic Net using stability Selection 
# ============================================================

# --- Setup ------------------------------------------------------------------
setwd()

# --- Load Required Packages -------------------------------------------------
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(paradox)
library(data.table)
library(mlr3misc)  # for some helper functions
library(mlr3mbo)
library(readr)
library(mlr3viz)
library(parallelMap)
library(parallel)
library(dplyr)
library(glmnet)

#load("Selection.RData")
load("Selection_PH.RData")


autoplot(rr)
set.seed(123)  

# 1. Score each outer fold
score_dt <- rr$score(msr("regr.mse"))

# 2. Identify the best-performing outer fold (i.e. the one with the lowest MSE)
best_fold_index <- which.min(score_dt$regr.mse)
cat("Best outer fold index:", best_fold_index, "\n")

# 3. Extract the tuning result (i.e. best hyperparameters) from that outer iteration
best_hyperpars <- rr$learners[[best_fold_index]]$tuning_result$x
cat("Best hyperparameters from the best outer fold:\n")
print(best_hyperpars)

#learner_elastic <- rr$learners[[best_fold_index]]$model
learner_elastic$param_set$values <- best_hyperpars[[1]]


# ============================================================
# Stability Selection
# ============================================================

# Define a new task for stability selection
task_stab <- TaskRegr$new(id = "gwas_rf", backend = gwas_data_rf, target = "PH")

# Number of subsampling runs for stability selection
n_runs <- 1000

# Initialize a named vector to store the count of selections for each feature
feature_names <- task_stab$feature_names  
selection_counts <- setNames(rep(0, length(feature_names)), feature_names)

# Stability selection: subsample data and record selection frequency

for (i in seq_len(n_runs)) {
  # Create a subsample: use 50% of the observations
  subsample_indices <- sample(seq_len(task_stab$nrow), size = floor(0.5 * task_stab$nrow))
  
  # Create a new task with the subsampled data
  subsample_data <- task_stab$data(rows = subsample_indices)
  task_sub <- TaskRegr$new(id = paste0("sub_", i), backend = subsample_data, target = "PH")
  
  # Train the elastic net model on the subsample
  learner_elastic$train(task_sub)
  
  # Extract the coefficient vector (s should be tuned already)
  coef_matrix <- coef(learner_elastic$model, s = learner_elastic$param_set$values$s)
  
  # Remove the intercept 
  coefs <- coef_matrix[-1, , drop = TRUE]
  
  # Identify features selected in this run (nonzero coefficients)
  selected_features <- names(coefs)[abs(coefs) > 0]
  
  # Update selection counts for features selected
  selection_counts[selected_features] <- selection_counts[selected_features] + 1
  
  cat("Run", i, "completed\n")
}

# Compute the selection frequency for each feature
selection_frequency <- selection_counts / n_runs

stability_results <- data.frame(SNP = names(selection_frequency),
                                SelectionFrequency = as.vector(selection_frequency)) %>%
  arrange(desc(SelectionFrequency))
print(head(stability_results, 10))

# Set a threshold for selection frequency at 0.95
threshold <- 0.95
stable_features <- stability_results$SNP[stability_results$SelectionFrequency >= threshold]
cat("Number of stable features (selection frequency >=", threshold, "):", length(stable_features), "\n")
print(stable_features)


# ============================================================
# Merge with Marker Metadata and Save Results
# ============================================================
#learner_rf$importance
save.image("mlr3_GLMnet_PH2.RData")
myGM <- read.csv("../../../imputation_test/new_data/myGM_PH.csv")
#myGM <- read.csv("../../imputation_test/Gapit_PH//myGM.csv")

importance_df2 <- merge(avg_importance, myGM, by = "SNP")
write.table(importance_df2, "importance_GLMnet_PH2.txt", row.names = FALSE)
write.table(emp_threshold_perm,"emp_threshold_GLMnet_PH2.txt",row.names = F)
stability_results<-merge(stability_results,myGM)
write.table(stability_results,"stability_selection_feautres_PH2.txt",row.names = F)


