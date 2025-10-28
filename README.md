# Applying traditional and machine learning-based GWAS approaches for marker-trait identification in wheat

This repository contains **R scripts** and **input data** used for conducting both **traditional GWAS** and **machine learning–based GWAS (ML-GWAS)** as described in the corresponding publication.  
The purpose of this repository is to enable **replication** of the analyses and facilitate adaptation of the presented workflow for related genomic studies.

This project compares traditional **Genome-Wide Association Studies (GWAS)** with **machine learning–based GWAS** approaches.  
While traditional GWAS identifies marker–trait associations using linear or mixed models, the ML-based GWAS framework leverages regression algorithms to predict traits and extract the most informative genomic features.

The **ML-GWAS** approach implemented here includes:
- **Trait prediction** using regression algorithms  
- **Feature extraction** to identify informative SNPs  
- **Empirical cutoff estimation** for feature importance  
- **Implementation within the `mlr3` framework**

- ### Traditional GWAS
- Conducted using standard R scripts for marker–trait association testing following the respective tools manual.
- Input files follow typical GWAS conventions (HapMap or numeric genotype matrices, and phenotype tables).

### ML-based GWAS
- Implemented within the [`mlr3`](https://mlr3.mlr-org.com/) ecosystem in R.  
- Regression algorithms applied:
  - **Random Forest (RF)**
  - **Extreme Gradient Boosting (XGBoost)**
  - **Two-Stage Algorithm Based on Least Angle Regression and Random Forest (TSLRF) ** following Sun, J., Wu, Q., Shen, D. et al. TSLRF: Two-Stage Algorithm Based on Least Angle Regression and Random Forest in genome-wide association studies. Sci Rep 9, 18034 (2019). https://doi.org/10.1038/s41598-019-54519-x
  - **Elastic Net (EN)**
- **Stability selection** used to identify robust SNP features for EN following Meinshausen, N., & Bühlmann, P. (2010). Stability selection. Journal of the Royal Statistical Society Series B: Statistical Methodology, Volume 72, Issue 4, September 2010, Pages 417–473, https://doi.org/10.1111/j.1467-9868.2010.00740.x

