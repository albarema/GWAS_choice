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
rule all:
    input:
        expand("", pheno=PHENOS) 

## global parameters


checkpoint polyAdapt_freqs:
    input:
        infile=config['ukb_file'],
        popfile=config['acf_file'],
        lbd=config['lbd_file']
    output:
        outmp=temp("data/{panel}/gwasfreqs-{level}-{pheno}.tsv"),
        outfile="data/{panel}/gwasfreqs-{level}-{pheno}.tsv.gz",
        candi="data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
        neut="data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv"
    params:
        pvalCan= 5e-8,
        pvalNeut= 1e-5
    shell:
        """
        python scripts/acf2ukbfreq_byP_Alba.py  -a {input.popfile} -g {input.infile} -o {output.outmp}
        cat <(head -1 {output.outmp}) <(tail -n+2 {output.outmp} | sort -k1,1 -k2,2g) | bgzip -c > {output.outfile}
        tabix -s 1 -b 2 -e 2 {output.outfile}
        python scripts/partitionUKB_byP.py -i {output.outfile} -b {input.lbd} -o {output.candi} -p {params.pvalCan}
        python scripts/extractneutral_byP.py -i {output.outfile} -o {output.neut} -s 20 -p {params.pvalNeut}
        """

rule polyAdapt_qx:
    input:
        neut="data/gwasfreqs_neutral-{level}-{pheno}.tsv",
        candi="data/gwasfreqs_candidates-{level}-{pheno}.tsv",
    output:
        qx="selection_gwas/QX_report-{level}-{pheno}.txt",
        scores="selection_gwas/Genscores-{level}-{pheno}.txt"
    shell:
        """
        Rscript scripts/CalcQX_edit4parallel_Alba.R -w {input.candi} -e {input.neut} -o {output.qx} -s {output.scores} -n 1000 -j 1000
        """
