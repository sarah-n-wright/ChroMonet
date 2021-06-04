# ChroMonet : *An HHM-based tool for local ancestry inference and chromosome painting*
Final class project for UCSD's CSE284 Personal Genomics, Spring 2021  
Authored by: Lauryn Bruce, Hannah Mummey & Sarah Wright

## Overview & Key Functionalities

**ChroMonet** is an Hidden Markov Model implementation for local ancestry inference. This implementation is designed to work with two ancestral populations, and SNP positions with base pair substitutions only (no indels). The model makes the simplifying assumption that recombination is equally likely at all positions in the genome.

This respository contains all necessary functions to implement, run and evaluate ChroMonet, as well as methods for creation of simulated datasets, genomes and admixtures. Note that the repository contains toy simulated datasets, but additional data needs to be created or downloaded by users to use ChroMonet with longer datasets.

## Data used to develop and test ChroMonet

### Simulated data

* `simulated_data.py` contains methods to generate simulated allele frequencies for a specified number of positions and produce synthetic genomes from these simulated datasets, and perform admixture between genomes. In all functions, the individual positions are considered independent, as are the two chromosomes within a genome. 

### 1000 Genomes Data

* `parse_genotype_mafs.py` Extraction and formatting of allele freqencies for AFR and EUR super-populations from 1000 Genomes VCF files. 

### Data Availability

`Data/` contains example simulated datasets for N=100 positions.   

Larger simulated datasets and 1000 Genomes data used in this project can be found on the class JupyterHub in Team3 directory under `simulated_files` and `chromosome_21_files`/`chromosome_14_files` respectively. Alternatively, simulated datasets can be generated as described below, and 1000 Genomes data can be download from [The 1000 Genomes Project](https://www.internationalgenome.org/).

## HMM implementation and evaluation

### Implementation
* `LAI_hmm_script.py` All methods for implementing local ancestry HMM and solving for optimal path. 
* `run_LAI_hmm.py` Example script based implementation of HMM.

### Evaluation
`accuracy_metrics.py` contains methods for evaluating the predicted results from the HMM by comparing the results to corresponding truth sets. Wrapper function `accuracy_metrics.compare_genome_files()` takes file paths to the predicted and truth genotypes and outputs:
  * Positional accuracy: Proportion of all alleles assigned the correct ancestry.
  * Skew: The ratio of number predicted recombinations events to the true number of recombination events.   

### Results Visualization

* `plot_karyogram_LB.py` Implementation of chromosome painting adopted from a script originally written by Alicia Martin  https://github.com/armartin/ancestry_pipeline. Modified to work with single chromosomes and simulated genomes. 
* `plot_chromosome_painting.py` Wrapper script for `plot_karyogram_LB.py` that formats HMM outputs for chromosome paintings and calls chromosome painting functions.

## Implementation Notebooks
These jupyter notebooks provide example implementation of the HMM and associated tasks. These are for illustration purposes, and will not run correctly outside of the CSE284 JupyterHub. 

### Simulated Data & HMM Usage
* `Simulated Data Examples.ipynb` Example usage of functions to create synthetic allele frequency datasets, synthetic genomes and perform admixture of synthetic and real genomes. 
* `File Reformatting.ipynb` Formalization of input data formats
* `Streamlined HMM.ipynb` Implementation and testing of HMM on simulated data
* `Log10 HMM.ipynb` Implementation and testing of HMM with model that utilizes addition of logs of probabilities rather than product of probabilities
* `HMM Performance Testing.ipynb` HMM implementation that outputs runtime estimates for runtime testing

### Evaluation
* `Accuracy calculations.ipynb` Calculation and visualization of accuracy metrics for all HMM tests
* `Create Chromosome Paintings.ipynb` Visulization of HMM outputs using chromosome paintings. 

### 1000 Genomes Data
* `Generate 1000G AFR-EUR Inputs.ipynb` Calculation of necessary HMM inputs from 1000 Genomes data
* `CHR14 AF Calculations.ipynb` Allele frequency calculations for Chromosome 21 in 1000 Genomes data
* `CHR21 AF Calculations.ipynb` Allele frequency calculations for Chromosome 14 in 1000 Genomes data


## File naming conventions

### Simualted data files
All simulated data files start with the prefix `sim`.   
* Simulated allele frequencies: `simData_N[number of populations]_P[number of positions]_seed[random seed].tsv`. E.g. the file `simData_N2_P100_seed518.tsv` is a two population data set with 100 positions, generated with random seed=518.
* Simulated genomes: `simGenome_[number of positions]_[pop1]_[pop2].tsv`. E.g. `simGenome_100_1_0.tsv` is a 100 position genome with one chromosome from population 1, and one chromsome from population 0.

### Admixed genomes
Simulated admixed genomes: `simAdmixedGenome_[number of positions]_[pop1]_[pop2]_Rx[number of recombinations]_[repeat].tsv` e.g. `simAdmixedGenome_100_0_1_Rx4_1.tsv` specifies the first repeat of admixture between population 0 and population 1 with 4 recombination events in each chromosome. 

### HMM output files
`[genome file prefix]_recomb[recombination rate]_HMMoutput.tsv` e.g. `simAdmixedGenome_100_0_1_Rx1_0_recomb0.01_HMMoutput.tsv` specifies the HMM prediction results for genome `simAdmixedGenome_100_0_1_Rx1_0.tsv` with recombination rate r=0.01
