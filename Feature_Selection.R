# ============================================================
# Feature selection using Elastic Net
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
#geno <- read.csv("../../../imputation_test/new_data/myGD_PH.hmp.txt")
geno <-read.csv("../../../imputation_test/new_data/myGD_TKW.hmp.txt")
#pheno <- read.table("../../../data/Cimmyt/PH_LODG.txt",header = T)
pheno <- read.csv("../../../data/Cimmyt/Qtraits.csv")
colnames(pheno)[1]<-"taxa"
# combine pheno & geno data
gwas_data <- merge(pheno[,c(1,5)], geno, by = "taxa")
#gwas_data <- merge(pheno[,c(1,3)], geno, by = "taxa")
gwas_data<-gwas_data[,-1]
# --- Regression Task --------------------------------------------------------

task_gwas = TaskRegr$new(id = "gwas", backend = gwas_data, target = "TKW")
#task_gwas = TaskRegr$new(id = "gwas", backend = gwas_data, target = "PH")
# ============================================================
# 1. Elastic Net Learner & Hyperparameter Tuning
# ============================================================
#as.data.table(mlr_tuning_spaces)
#tuning_space = lts("regr.glmnet.default")
#tuning_space
#as.data.table(lrn("regr.glmnet")$param_set)

learner_elastic = lrn("regr.glmnet",parallel=TRUE,
                      s = to_tune(1e-5, 10000, logscale =T),# penalize complex models
                      alpha = to_tune(0, 1))# decicide if ridge regression, lasso or mix 

# Define tuner, measure, and terminator
tuner = tnr("mbo") # bayesian optemization
#tuner = tnr("random_search")
measure_mse = msr("regr.mse")
terminator_eval =trm("evals", n_evals = 10)

# Define resampling strategies
outer_rs <- rsmp("cv", folds = 5)  # outer loop
inner_rs <- rsmp("cv", folds = 5)  # inner loop

# Nested CV tuner
at = auto_tuner(
  tuner = tuner,
  learner = learner_elastic,
  resampling = inner_rs,
  measure = msr("regr.mse"),
  terminator = terminator_eval
)
# --- Run Nested Cross-Validation --------------------------------------------
rr = mlr3::resample(task_gwas, at, outer_rs,store_models = T)

# -- Select Features ----------------------------------------------------------
coef_list <- list()
for (model in 1:5) {
  mod<-rr$learners[[model]]$model
  coefs = as.vector(coef(mod$learner$model, s = mod$tuning_instance$result_learner_param_vals$s))
  names(coefs) = rownames(coef(mod$learner$model, s = mod$tuning_instance$result_learner_param_vals$s))
  coefs = coefs[-which(names(coefs) == "(Intercept)")]
  coef_list[[model]] <- coefs
}

coef_dt <- rbindlist(lapply(coef_list, function(x) as.data.table(as.list(x))), fill = TRUE)

# Replace NA values with 0 (for features not selected in some models)
coef_dt[is.na(coef_dt)] <- 0

# Compute the mean coefficient for each feature
mean_coefs <- coef_dt[, lapply(.SD, mean)]

# Transpose the data.table 
mean_coefs<-as.data.frame(t(mean_coefs))

# Filter out features with a mean coefficient of zero
selected_features <- which(mean_coefs$V1!= 0)
mean_coefs2<-mean_coefs[selected_features,,drop=F]
selected_features2<-rownames(mean_coefs2)

# ============================================================
# 2. Save data
# ============================================================

selected_cols = c("PH", selected_features2)
gwas_data_rf = gwas_data[, selected_cols, drop = FALSE]

write.table(gwas_data_rf,"Selected_SNP_PH.txt",col.names = T,row.names = F,quote = F)
save.image("Selection_PH.RData")
