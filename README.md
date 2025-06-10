# Forecasting Count Data with Varying Dispersion: A Latent-Variable Approach

This repository contains code for reproducing the simulation study and case studies for the forthcoming paper "Forecasting Count Data with Varying Dispersion: A Latent-Variable Approach" by Easton K. Huch, Candace Berrett, Mason Ferlic, and Kimberly F. Sellers.

To reproduce the figure comparing the approximate expected Fisher information for the GP-P model to its exact value, execute the script `eim_numerical_integration.R`.

To reproduce the prediction interval simulation, execute the script `pred-interval-simulation.R`.

To reproduce the main simulation study, execute the R files in the `simulation_study` directory.

To reproduce the Dominicks case study, first download the data files for the `cer`, `coo`, and `cra` categories into the corresponding directories within the `dominicks-data` directory.
This can be accomplished by rendering the notebook `download-dominicks-data.ipynb` or by manually downloading the files from the [Kilts Center for Marketing](https://www.chicagobooth.edu/research/kilts/research-data/dominicks).
Second, run the script `create-dominicks-data-set.R` to perform preprocessing.
Third, run the script `dominicks.R` to analyze the processed data.

To reproduce the COVID case study, download the files `country-lookup.csv`, `index.csv`, `demographics.csv`, `economy.csv`, `health.csv`, and `epidemiology.csv` from [Google Health](https://health.google.com/covid-19/open-data/raw-data).
Then execute the file `create-covid-data-set.R`.
Optionally, execute the file `covid-model-selection.R` to reproduce the model selection process.
Finally, execute the file `covid-case-study.R`.