```diff
- NB! CURRENTLY UNDER DEVELOPMENT 
```

# Evaluating the robustness of polygenic adaptation to the choice of the GWAS cohort


# A - Neutrality test for polygenic scores

Summary files of 30 phenotypes shared at least between 2 GWAS: 
- phenotypes.both_sexes.tsv.gz

Genomic data:
- Run vcf2acf.smk

Compute Qx statistic
- Run polyadapt.smk 

# B - Assessing different association methods
Quality control filters:
- Genotyped SNPs
- MAF > 0.001
- HWE pvalue > 1e-10
- Autosomes

Association models - examples

- Linear model implemented in PLINK 2.0
```bash 
plink2
```
- Linear mixed model implemented in BOLT_LMM
```bash 
plink2
```
- Meta-analysis implemented in METAL
```bash 
//inverse variance based
metal

// sample size based
metal
```


