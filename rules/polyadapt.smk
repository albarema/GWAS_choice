#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:01:16 2020
@author: gsd818
"""

import pandas as pd
import csv

configfile: "config.yaml"

PHENOS=pd.read_table('phenotypes.filtered.both_sexes.txt')['phenotype'].tolist()
## rules
rule betas:
    input:
        expand(os.path.join(config['uk_dir'], "{pheno}.flipped.byP.gz"), pheno=PHENOS) 

## global parameters


checkpoint polyAdapt_freqs:
    input:
        infile=config['ukb_file'],
        popfile=config['acf_file'],
        lbd=config['lbd_file']
    output:
        outfile="data/{panel}/gwasfreqs-{level}-{pheno}.tsv.gz",
        candi="data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv"
    shell:
        """
        python scripts/acf2ukbfreq_byP_Alba.py  -a {input.popfile} -g {input.infile} -o UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv
        cat <(head -1 UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv) <(tail -n+2 UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv| sort -k1,1 -k2,2g) | bgzip -c > {output.outfile}
        tabix -s 1 -b 2 -e 2 {output.outfile}
        python scripts/partitionUKB_byP.py -i {output.outfile} -b {input.lbd} -o {output.candi} -p 5e-08
        python scripts/extractneutral_byP.py -i {output.outfile} -o {output.neut} -s 20 -p 0.00001
        """

rule polyAdapt_qx:
    input:
        neut="UKBiobank/data/gwasfreqs_neutral-{level}-{pheno}.tsv",
        candi="UKBiobank/data/gwasfreqs_candidates-{level}-{pheno}.tsv",
    output:
        qx="UKBiobank/selection_UKBV2/QX_report-{level}-{pheno}.txt",
        scores="UKBiobank/selection_UKBV2/Genscores-{level}-{pheno}.txt"
    shell:
        """
        Rscript scripts/CalcQX_edit4parallel_Alba.R -w {input.candi} -e {input.neut} -o {output.qx} -s {output.scores} -n 1000 -j 1000
        """

rule qx_distr:
    input:
        "phenoname_qx.txt"
    output:
        allqx="UKBiobank/Qx-pvalue/Qx_allvals_alltraits.txt",
        sigqx="UKBiobank/Qx-pvalue/Qx_allvals_sigtraits.txt"
    shell:
        """
        Rscript scripts/QX_distr_plots.R {input} {output.allqx} {output.sigqx}
        """


rule pol_scores:
    input:
        cols="clusters_name_cols.txt",
        qx="UKBiobank/Qx-pvalue/Qx_allvals_sigtraits.txt",
        categ="traits_descript_categories.tsv"
    output:
        "UKBiobank/polyscores_sigtraits_5e-8.txt"
    shell:
        """
        Rscript scripts/scoresPlot.R {input.qx} {input.cols} {input.categ}
        """






