# Romanelli et al. 2021

This repo contains the data and analysis scripts used to generate the RNA-seq figures in Romanelli et al. *J. Biol. Chem* 2021 (full citation below). The raw FASTQ files and gene counts from `STAR` can be found at GEO ([GSE176453](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176453)) or in the **alignments** directory here. The **scripts** directory contains additional mapping scripts that are for reference only and not required to replicate the results.

## Files
```
+-- alignments
    +-- 1131BAT
        +-- ReadsPerGene.out.tab
    +-- 1177BAT
        +-- ReadsPerGene.out.tab
    +-- ...
+-- data
    +-- genes.csv
    +-- samples.csv
+-- scripts
    +-- analysis.R
+-- environment.yml
```

## Pipeline
1. Create conda environment: `conda env create -f environment.yml`
2. Activate `conda` environment: `conda activate romanelli`
3. Run the analysis and generate outputs: `Rscript scripts/analysis.R`
4. Deactivate the `conda` environment: `conda deactivate`

## Citation
(full citation to come)
