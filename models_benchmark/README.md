# TCR-pMHC Specificity Benchmark

This repository contains notebooks used to perform the calculations for the ML TCR-pMHC specificity benchmark as described in the paper. Four models (ATM-TCR, Gradient Boosting, Random Forest and TULIP) are benchmarked using the dataset of TULIP (Meynard-Piganeau et al., 2024). The code of Random Forest and Gradient Boosting is available in this repository, ATM-TCR and TULIP can be run from their public release.



## Environment Setup

To set up the environment for the Random Forest and XGBoost models, use the provided environment yaml file. Follow these steps:

1. Ensure you have Conda installed on your system. If not, you can download and install it from [here](https://docs.conda.io/en/latest/miniconda.html).
2. Clone this repository:
    ```bash
    git clone https://github.com/instadeepai/tcr-dataset-size-estimation/
    cd tcr-dataset-size-estimation/models_benchmark/
    ```
3. Create the Conda environment from the YAML file:
    ```bash
    conda env create -f environment.yml
    ```
4. Install the package in editable mode (to avoid issues with Hydra):
    ```bash
    pip install -e .
    ```
5. Activate the environment:
    ```bash
    conda activate tcr_pmhc_benchmark_env
    ```


## Description of the notebooks

### 1. Data Preparation (`dataset_preparation.ipynb`)
- Loads the data from sourced from Meynard-Piganeau et al., 2024 (TULIP model).
- Applies several filters to the dataset to remove:
  1. Duplicates.
  2. Non-MHC-I samples.
  3. Missing sequences.
  4. Too long CDR3b / peptide sequences.
  5. Samples for which a CDR3b is paired to a train peptide and a test peptide, both being at edit distance 1.
  6. Peptides with less (resp. more) than 20 samples in the test (resp. train) set.

### 2. Nested Bootstrap and ROC-AUC Calculation (`bootstrap.ipynb`)
- Collects the results of the models.
- Performs a nested bootstrap to obtain the ROC-AUC confidence intervals of each model as a function of the peptide edit distance to the train set.
- To run this notebook for each model benchmarked, change the variable `model` at the beginning of the notebook.

### 3. ROC-AUC/PR-AUC barplots (`plot_fig2a_2b.ipynb`)
- Draws the barplot for Figure 2a 2b in the paper.
- Represents the ROC-AUC and PR-AUC per model.

### 4. ROC-AUC per peptide edit distance (`plot_fig2c.ipynb`)
- Represents the ROC-AUC as a function of the peptide edit distance to the train set for each model using the confidence intervals obtained in Notebook N2.
- Corresponds to Figure 2c in the paper.

## Usage

To use the code and notebooks, follow these steps (you may skip to steps 4-5 using the data provided):
1. **We do not distribute public databases with this work. TULIP dataset should be accessed from the [TULIP GitHub repository](https://github.com/barthelemymp/TULIP-TCR/blob/main/data/).**
2. Open the notebook `dataset_preparation.ipynb` to prepare the data.
3. Run the models:  
    - **Random Forest**  
      Generate the embeddings  
    ```bash
    embedding_generation data.csv_path=dataset/models_benchmark_dataset.csv data.data_dir=rf_model_benchmark
    ```
      Train and run the Random Forest  
    ```bash
    train_predict data.data_dir=rf_model_benchmark data.output_dir=results_model model=rf
    ```
   - **Gradient Boosting**  
    Generate the embeddings  
    ```bash
    embedding_generation data.csv_path=dataset/models_benchmark_dataset.csv data.data_dir=xgb_model_benchmark
    ```
      Train and run the Gradient Boosting  
    ```bash
    train_predict data.data_dir=xgb_model_benchmark data.output_dir=results_model model=xgb
    ```
	- **TULIP**  
  Clone [TULIP GitHub repository](https://github.com/barthelemymp/TULIP-TCR/).  
  Instructions to run the TULIP model are provided in the TULIP README file. Input and output dataframes need to be reformatted appropriately.
	- **ATM-TCR**  
  Clone [ATM-TCR GitHub repository](https://github.com/Lee-CBG/ATM-TCR).  
  Copy the `atm_tcr.patch` file located in `models_benchmark/models` into your local cloned ATM-TCR repository and apply the patch to add the focal loss.
    ```bash
    patch -p1 <atm_tcr.patch
    ```
      Follow the instructions provided in the ATM-TCR README file to run it.  Input and output dataframes need to be reformatted appropriately.

4. Open the notebook `bootstrap.ipynb`, set the `model` variable to the desired model, and run the notebook to perform the nested bootstrap analysis.
5. Open the notebook `plot_fig2a_2b.ipynb` to generate the barplot for the benchmark results.
6. Open the notebook `plot_fig2c.ipynb` to visualize the ROC-AUC as a function of the peptide edit distance to the train set for each model.


## References
```
@article{TULIP,
	author = {Meynard-Piganeau, Barthelemy and Feinauer, Christoph and Weigt, Martin and Walczak, Aleksandra M. and Mora, Thierry},
	title = {TULIP {\textemdash} a Transformer based Unsupervised Language model for Interacting Peptides and T-cell receptors that generalizes to unseen epitopes},
	year = {2024},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/02/02/2023.07.19.549669},
	journal = {bioRxiv}
}
```
```
@article{ATM-TCR,
  title={ATM-TCR: TCR-epitope binding affinity prediction using a multi-head self-attention model},
  author={Cai, Michael and Bang, Seojin and Zhang, Pengfei and Lee, Heewook},
  journal={Frontiers in immunology},
  volume={13},
  year={2022},
  publisher={Frontiers}
}
```
## License
- The code is distributed under the GPL v3 license, accessible [here](LICENSE)  
- ATM-TCR is distributed under CC-BY-4.0 license, accessible [here](https://github.com/Lee-CBG/ATM-TCR?tab=CC-BY-4.0-1-ov-file).  
- TULIP is distributed under GPL-3.0 license, accessible [here](https://github.com/barthelemymp/TULIP-TCR?tab=GPL-3.0-1-ov-file).  
- Scikit-learn is distributed under BSD-3-Clause 3.0 license, accessible [here](https://github.com/scikit-learn/scikit-learn?tab=BSD-3-Clause-1-ov-file).  

## Citation
Please refer to the following:

```
@article{Data set size requirements for generalizable TCR-antigen specificity prediction,
  title={Assessing data size requirements for training generalizable sequence-based models for TCR-antigen specificity based on pan-allelic MHC-I non-self ligandome estimation},
  author={Delaunay, A. and McGibbon, M. and Djermani, B. and Gorbushin, N. and Chaves García-Mascaraque, S. and Rayment, I. and Kizhvatov, I. and Petit, C. and Lang, M. and Rooney, M. and Beguir, K. and Sahin, U. and Mikhaylov, V. and Copoiu, L. and Lopez Carranza, N. and Tovchigrechko, A.},
  journal={},
  volume={},
  year={2024},
  publisher={}
}
```

