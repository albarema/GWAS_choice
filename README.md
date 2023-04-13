```diff
- NB! CURRENTLY UNDER DEVELOPMENT 
```

# Evaluating the robustness of polygenic adaptation to the choice of the GWAS cohort


# A - Neutrality test for polygenic scores 
> Workflow 
[![INSERT YOUR GRAPHIC HERE](workflow.example.png)]()

> Step 1: download all files 
You can modify the config.yaml with your own paths and dirs
For GRCh37, download files from:
- epo_file: 
```bash 
wget ftp://ftp.healthtech.dtu.dk/public/EPO/all_hg19.epo.gz; wget ftp://ftp.healthtech.dtu.dk/public/EPO/all_hg19.epo.gz.tbi
```
- fai_file: 
```bash 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
```
- lbd_blocks:
```bash 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
```
> Step 2: Genomic data - get your vcf ready
You can get the allele frequency of the populations of interest here (e.g.: populations from the 1000 Genomes Project) by running vcf2acf.smk using snakemake 

```bash 
snakemake --snakefile path/to/vcf2acf.smk --cores XX
```

More info on how to use glactools here: https://github.com/grenaud/glactools

> Data 
- Summary files of 30 phenotypes shared at least between 2 GWAS


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
plink2 --bgen allind.bgen --sample allind.sample --hwe 1e-10 --maf 0.001 --keep-fam filtered_ind.list.txt 
--linear hide-covar --pheno $HEIGHT --covar allcov.tsv --covar-col-nums 3-4,8-27 --variance-standardize --out $WHITE
```
- Meta-analysis implemented in METAL
```bash 
//inverse variance based
metal

// sample size based
metal
```


