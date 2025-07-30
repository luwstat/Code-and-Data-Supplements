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
