# Applying traditional and machine learning-based GWAS approaches for marker-trait identification in wheat

This repository contains **R scripts** used for conducting both **traditional GWAS** and **machine learning–based GWAS (ML-GWAS)** as described in the corresponding publication.  
The purpose of this repository is to enable **replication** of the analyses and facilitate adaptation of the presented workflow for related genomic studies.

This project compares traditional **Genome-Wide Association Studies (GWAS)** with **machine learning–based GWAS** approaches.  
While traditional GWAS identifies marker–trait associations using linear or mixed models, the ML-based GWAS framework leverages regression algorithms to predict traits and extract the most informative genomic features.

The **ML-GWAS** approach implemented here includes:
- **Trait prediction** using regression algorithms  
- **Feature extraction** to identify informative SNPs  
- **Empirical cutoff estimation** for feature importance  
- **Implementation within the [`mlr3`](https://mlr3.mlr-org.com/) framework**

### Data source
The data analyzed in this study were obtained from:

> **Juliana, P., Poland, J., Huerta-Espino, J., Shrestha, S., Crossa, J., Crespo-Herrera, L., Toledo, F. H., Govindan, V., Mondal, S., Kumar, U., Bhavani, S., Singh, P. K., Randhawa, M. S., He, X., Guzman, C., Dreisigacker, S., Rouse, M. N., Jin, Y., Pérez-Rodríguez, P., … Singh, R. P. (2019).**  
> *Improving grain yield, stress resilience and quality of bread wheat using large-scale genomics.*  
> **Nature Genetics, 51(10), 1530–1539.**  
> [https://doi.org/10.1038/s41588-019-0496-6](https://doi.org/10.1038/s41588-019-0496-6)

The dataset includes genotypic and phenotypic information from a large, diverse panel of bread wheat accessions and was used as the foundation for traditional and ML-based GWAS analyses in this project.

### Input Files

- **Traditional GWAS ([`SOMMER`](https://www.rdocumentation.org/packages/sommer/versions/4.4.3) and [`GAPIT`](https://zzlab.net/GAPIT/)):**  
  Input files should follow the formats described in the respective package manuals.  
  Ensure that genotype and phenotype data are properly aligned by sample identifiers.

- **ML-based GWAS:**  
  The input file should begin with the **trait of interest** (phenotypic values), followed by the **marker data** in **numeric format**.  
  - Marker names should be used as **column headers**.  
  - Individual (sample) names should be used as **row identifiers**.  
  - Do **not** include chromosome or position information in this file — only marker values.  
  - Each cell should represent the **numerical genotype code** (e.g., 0, 1, 2) or (-1, 0, 1).


- ### Traditional GWAS
- Conducted using standard R scripts for marker–trait association testing following the respective tools manual.

### ML-based GWAS
- Implemented within the [`mlr3`](https://mlr3.mlr-org.com/) ecosystem in R.  
- Regression algorithms applied:
  - **Random Forest (RF)**
  - **Extreme Gradient Boosting (XGBoost)**
  - **Two-Stage Algorithm Based on Least Angle Regression and Random Forest (TSLRF)** following Sun, J., Wu, Q., Shen, D. et al. TSLRF: Two-Stage Algorithm Based on Least Angle Regression and Random Forest in genome-wide association studies. Sci Rep 9, 18034 (2019). https://doi.org/10.1038/s41598-019-54519-x
  - **Elastic Net (EN)**
- **Stability selection** used to identify robust SNP features for EN following Meinshausen, N., & Bühlmann, P. (2010). Stability selection. Journal of the Royal Statistical Society Series B: Statistical Methodology, Volume 72, Issue 4, September 2010, Pages 417–473, https://doi.org/10.1111/j.1467-9868.2010.00740.x

