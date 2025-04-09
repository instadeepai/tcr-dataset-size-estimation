# Assessing Data Size Requirements For Training Generalizable Sequence-based Models of TCR-antigen Specificity Based on Pan-allelic MHC-I Non-self Ligandome Estimation

This repository contains the companion code for the paper [Assessing data size requirements for training generalizable sequence-based models for TCR-antigen specificity based on pan-allelic MHC-I non-self ligandome estimation]() by Delaunay A., McGibbon M et al.

# Usage

The repository follows the structure of the paper:

* The `models_benchmark` folder contains all code relevant for preparing the benchmark dataset, calculating bootstrapped ROC-AUC of models' specificity prediction and reproducing the paper's visualization plots.

* The `ligandome` folder contains the code for reproducing the ligandome computation workflow on public data. 

Please refer to each folder's README file for environment setup and specific usage instructions.

## License and Copyright

This code is copyright of BioNTech SE, 2022-2024.

The code is available under the [GPL v3](LICENSE) license terms.
Due to licensing restrictions regarding redistribution for some of the external datasets, we provide instructions within this repository to access and prepare all datasets used in this work from their original sources.

## Citing this work

Please refer to the following:

```
@article{Data set size requirements for generalizable TCR-antigen specificity prediction,
  title={Assessing data size requirements for training generalizable sequence-based models for TCR-antigen specificity based on pan-allelic MHC-I non-self ligandome estimation},
  author={Delaunay, A. and McGibbon, M. and Djermani, B. and Gorbushin, N. and Chaves García-Mascaraque, S. and Rayment, I. and Kizhvatov, I. and Petit, C. and Lang, M. and Rooney, M. and Beguir, K. and Sahin, U. and Copoiu, L. and Lopez Carranza, N. and Tovchigrechko, A.},
  journal={},
  volume={},
  year={2024},
  publisher={}
}
```
