# Code-and-Data-Supplements
This repository contains all code and simulated data necessary to reproduce the simulation results presented in our paper titled "Regression Analysis of Arbitrarily Censored and Left-Truncated Data Under the Proportional Odds Model." Please refer to the README.txt file to know more details.


### 1. Simulated Datasets

The `Simulated_Datasets` folder contains all datasets generated and used in the simulation section of the manuscript.

- Each `.txt` file represents a dataset with 500 replicates, each of size 200.
- The filenames encode key information:
  - The **last two numbers** in the filename indicate the **true values of the regression coefficients** β.
  - If the filename contains the letter **"r"**, the data were generated under **arbitrary censoring without left truncation**.
  - If **"r"** is not present, the dataset corresponds to **arbitrarily censored and left-truncated** data.

These datasets are used to evaluate and compare model performance under different censoring and truncation scenarios. Each dataset includes the necessary variables for fitting the proportional odds model and reproducing the simulation results reported in the manuscript.


### 2. Simulation Code

The `Simulation_Code` folder contains all the scripts used to fit models and summarize results automatically.

- `em_fit.R`: Fits the model using the **proposed EM algorithm**.
- `icen_fit.R`: Fits the model using the `ic_sp` function from the **icenReg** package.
- `project1EMfit.R`: Fits the model using the **comparison method** implemented in the `regPO` package.
- `do_rt_auto_icen.R`: Automatically summarizes results from `icen_fit.R` and generates the output tables used in the manuscript.
- `do_rt_auto_trun.R`: Automatically summarizes results from `em_fit.R` and produces the corresponding tables for the manuscript.
- `do_rt_auto.R`: Automatically summarizes results from `project1EMfit.R` and generates the manuscript-ready tables.

These scripts are organized to enable full reproducibility of the simulation section with minimal manual intervention. Output files are saved with descriptive filenames and match the tables and figures presented in the manuscript.


### 3. Real Data Analysis Code

This folder contains the R code used to produce the results presented in the **real data analysis** section of the manuscript.

> **Note:** Due to data sharing restrictions, we are unable to include the raw datasets. Interested readers may request access to the data from the authors of the original studies:

1. Peeters, M., Price, T. J., Cervantes, A., Sobrero, A. F., Ducreux, M., Hotko, Y., André, T., Chan, E., Lordick, F., Punt, C. J. A., Strickland, A. H., Wilson, G., Ciuleanu, T.-E., Roman, L., Van Cutsem, E., Tzekova, V., Collins, S., Oliner, K. S., Rong, A., Gansert, J., et al. (2010). *Randomized Phase III study of panitumumab with fluorouracil, leucovorin, and irinotecan (FOLFIRI) compared with FOLFIRI alone as second-line treatment in patients with metastatic colorectal cancer.* Journal of Clinical Oncology, 28(31), 4706–4713. https://doi.org/10.1200/JCO.2009.27.6055

2. Tu, W., Batteiger, B. E., Wiehe, S. E., Ofner, S., Van Der Pol, B., Katz, B. P., Orr, D. P., & Fortenberry, J. D. (2009). *Time from first intercourse to first sexually transmitted infection diagnosis among adolescent women.* Archives of Pediatrics & Adolescent Medicine, 163(12), 1106–1111. https://doi.org/10.1001/archpediatrics.2009.203


#### **First Real Data Application (Metastatic Colorectal Cancer)**

- `colorectal_split.R`: Splits the dataset based on **KRAS mutation status** and fits the model using the proposed method.
- `colorectal_time_visual.R`: Visualizes the **time intervals** and **censoring types** of the survival outcome.
- `mCRC_auto_selection.R`: Runs model selection by iterating over different numbers of **interior knots** and **spline degrees**, and returns the corresponding **AIC** and **BIC** values.

#### **Second Real Data Application (Sexually Transmitted Infections)**

- `STD_TV.R`: Fits the model using the proposed method.
- `std_auto_degree_knots.R`: Performs model selection by varying the number of **interior knots** and **spline degrees**, and returns **AIC** and **BIC** values.
- `visual_noofpartner.R`: Visualizes the distribution of one covariate: **number of lifetime partners at enrollment**.
- `visual_survival_intervals.R`: Visualizes the **left-truncated survival intervals** and **censoring types** of the response variable.

Each script is clearly labeled and organized for easy execution and interpretation of the corresponding results in the manuscript.

