```diff
- NB! CURRENTLY UNDER DEVELOPMENT 
```

# Evaluating the robustness of polygenic adaptation to the choice of the GWAS cohort


# A - Neutrality test for polygenic scores
> Workflow 
[![INSERT YOUR GRAPHIC HERE](workflow.example.png)]()

> Data 
- Summary files of 30 phenotypes shared at least between 2 GWAS
- Genomic data: 
  - Allele frequencies from the 1000 Genomes Project
    - Run vcf2acf.smk

> Compute Qx statistic
- Before computing polygenic scores, if you are working with more than one GWAS, make sure their effect size is polarized by the same allele (i.e.: dervied allele)
  -  Run polyadapt.smk 

# B - Assessing different association methods
> Quality control filters: 

- Genotyped SNPs (~ 800,000)
- MAF > 0.001
- HWE pvalue > 1e-10
- Autosomes

> Association models
**Examples**
- Linear model implemented in PLINK 2.0
```bash 
plink2
```
- Linear mixed model implemented in BOLT_LMM
```bash 
WHITE="/home/white"
BRI="/home/bri"
HEIGHT="/home/id_height.tsv" # tabulated file with #FID	IID	HEIGHT

cd /home/white
plink2 --bgen allind.bgen --sample allind.sample --hwe 1e-10 --maf 0.001 --keep-fam filtered_ind.list.txt --linear hide-covar 
--pheno $HEIGHT --covar allcov.tsv --covar-col-nums 3-4,8-27 --variance-standardize --out $WHITE
```
- Meta-analysis implemented in METAL
```bash 
//inverse variance based
metal

// sample size based
metal
```


