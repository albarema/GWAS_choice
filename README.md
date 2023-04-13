```diff
- NB! CURRENTLY UNDER DEVELOPMENT 
```

# Evaluating the robustness of polygenic adaptation to the choice of the GWAS cohort


## A - Neutrality test for polygenic scores 
### Workflow 
[![INSERT YOUR GRAPHIC HERE](workflow.example.png)]()

### Step 0: download softwares and packages 
This workflow is built in snakemake. Please, download snakemake following the instructions here:
``` https://snakemake.readthedocs.io/en/stable/getting_started/installation.html ```

Install glactools to calculate alelle frequencies for a a given population panel. Follow instructions here:
``` https://github.com/grenaud/glactools ``` 

### Step 1: download all files 
You should modify the config.yaml with your own paths and dirs. If you are working with GRCh37, download files from:
 
```bash 
# EPO file
wget ftp://ftp.healthtech.dtu.dk/public/EPO/all_hg19.epo.gz; wget ftp://ftp.healthtech.dtu.dk/public/EPO/all_hg19.epo.gz.tbi

# FAI file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai

# LBD_blocks
# Download accordingly based on ancestry of the populations you are working with, e.g. EUR
wget https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/fourier_ls-all.bed

# Summary stats
# You can find a list of traits available for the UKBB by Neale lab at: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=178908679. For instance, for height: 
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O 50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz

# VCF file - 1000 Genomes can be download from here: 
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ 
```
Most files for GRCh38 are also available. However, you will need to liftover the LBD_blocks files. 

### Step 2: Genomic data and allele frequencies 
Get the allele frequency of the populations of interest (e.g.: populations from the 1000 Genomes Project) by running vcf2acf.smk using snakemake. The first rule is intended to get your vcf file ready (filters at the variant- and/or individual-level).

In the config file, add the path/filename under 'pops_panel'. This is a two-columnns tab-separate file with no header, first-column correspond to the sample ID and the second to the population they belong two. Example:
> head -n1 pops_panel_example_1000GP.txt
> 
> HG00096	GBR 

Once, the config file is ready. Run snakemake as follows:
```bash 
snakemake --snakefile path/to/vcf2acf.smk --cores XX
```

More info on how to use glactools here: https://github.com/grenaud/glactools

### Compute Qx statistic
Before computing polygenic scores, if you are working with more than one GWAS, make sure their effect size is polarized by the same allele (e.g.: dervied allele)
  -  Run polyadapt.smk 

## B - Assessing different association methods
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


