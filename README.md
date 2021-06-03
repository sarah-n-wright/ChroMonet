# ChroMonet : An HHM-based tool for local ancestry inference and chromosome painting
Final class project for UCSD's CSE284 Personal Genomics, Spring 2021
Authored by: Lauryn Bruce, Hannah Mummey & Sarah Wright

## Overview & Key Functionalities



## Data used to develop and test ChroMonet

### Simulated data

* `simulated_data.py` contains methods to generate simulated allele frequencies for a specified number of positions and produce synthetic genomes from these simulated datasets, and perform admixture between genomes. In all functions, the individual positions are considered independent, as are the two chromosomes within a genome. 

### 1000 Genomes Data

* `parse_genotype_mafs.py`

### Data Availability

`Data/` contains example simulated datasets for N=100 positions.   

Larger simulated datasets and 1000 Genomes data used in this project can be found on the class JupyterHub in Team3 directory under `simulated_files` and `chromosome_21_files`/`chromosome_14_files` respectively. Alternatively, simulated datasets can be generated as described below, and 1000 Genomes data can be download from [The 1000 Genomes Project](https://www.internationalgenome.org/).


## HMM implementation and evaluation

* `LAI_hmm_script.py`
* `run_LAI_hmm.py`

### Evaluation
* `accuracy_metrics.py` contains methods for evaluating the predicted results from and HMM by comparing the results to corresponding truth sets. Wrapper function `compare_genome_files()` takes file paths to the predicted and truth genotypes and outputs:
  * Positional accuracy: Proportion of all alleles assigned the correct ancestry.
  * Skew: The ratio of number predicted recombinations events to the true number of recombination events.   

### Results Visualization

* `plot_karyogram_LB.py`
* `plot_chromosome_painting.py`

## Interactive Notebooks
* `Simulated Data Examples.ipynb`
* `Create Chromosome Paintings.ipynb` 
* `generate_1000G_afr_eur_input_files.ipynb`


## File naming conventions

### Simualted data files
All simulated data files start with the prefix `sim`.   
* Simulated allele frequencies: `simData_N[number of populations]_P[number of positions]_seed[random seed].tsv`. E.g. the file `simData_N2_P100_seed518.tsv` is a two population data set with 100 positions, generated with random seed=518.
* Simulated genomes: `simGenome_[number of positions]_[pop1]_[pop2].tsv`. E.g. `simGenome_100_1_0.tsv` is a 100 position genome with one chromosome from population 1, and one chromsome from population 0.
* Simulated admixed genomes: `simAdmixedGenome_[number of positions]_[pop1]_[pop2]_Rx[number of recombinations]_[repeat].tsv` e.g. `simAdmixedGenome_100_0_1_Rx4_1.tsv` specifies the first repeat of admixture 

### Admixed genomes

### HMM output files

### Results visualization


