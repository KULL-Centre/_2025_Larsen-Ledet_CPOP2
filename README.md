# Blocking protein quality degradation leads to structural stabilization of DHFR indel variants
## Introduction
This respository contains all data (except from the raw FASTQ files, which are available at the NCBI Gene Expression Omnibus (GEO) repository (accession number: GSE302176)) and code to repeat the processing and analysis of the CPOP2 data in Larsen-Ledet et al.: "Blocking protein quality degradation leads to structural stabilization of DHFR indel variants".

## Overview of files
*Output files*
* **CPOP_data_WT.csv** - CPOP scores and standard deviations for DHFR indel, synonymous and nonsense variants in the WT strain.
* **CPOP_data_ubr1.csv** - CPOP scores and standard deviations for DHFR indel, synonymous and nonsense variants in the ubr1 knockout strain.
* **CPOP_data_san.csv** - CPOP scores and standard deviations for DHFR indel, synonymous and nonsense variants in the san1 knockout strain.
* **[WT|Ubr1|San1]_per_tile_variant_counts_tile[1-5].csv** - Counts per tile for DHFR indel, synonymous and nonsense variants for each replicate, condition, and strain (WT, ubr1 knockout, and san1 knockout).
  
*Input files*
* **CPOP_data_original_2024.csv** - CPOP scores and standard deviations for DHFR indel, synonymous and nonsense variants in the original CPOP paper Larsen-Ledet et al.: "Systematic characterization of indel variants using a yeast-based protein folding sensor" (https://doi.org/10.1016/j.str.2024.11.017)

*Excel files*
* **SupplementaryFile1.xlsx** - CPOP scores and standard deviations for DHFR indel, synonymous and nonsense variants in the WT strain, ubr1 knockout strain, and san1 knockout strain, as well as all primer sequences combined in a single Excel file.

## Processing of raw sequencing data
The function.py file (available here https://github.com/KULL-Centre/_2024_Larsen-Ledet_CPOP) is used to call DHFR variants and calculate CPOP scores. The script takes raw FASTQ files as input. The output is a dataset with CPOP scores and standard deviations for DHFR indel, synonymous and nonsense variants.

## Data analysis and plotting
The CPOP2.0_code.R file is used to produce all plots in the main figures, and the CPOP2.0_supplementary_code.R file is used to produce all plots in the supplementary figures.

## Paper

