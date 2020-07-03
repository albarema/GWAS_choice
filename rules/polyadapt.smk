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

rule polar_beta:
    input:
        os.path.join(config['uk_dir'],"{pheno}.gwas.imputed_v3.both_sexes.tsv.gz")
    output:
        ungz=temp(os.path.join(config['uk_dir'], "{pheno}.flipped.byP")),
        gz=os.path.join(config['uk_dir'], "{pheno}.flipped.byP.gz")
    shell:
        """
        python scripts/FlipMinorUKB.py -i {input} -o {output.ungz} 
        cat <(head -1 {output.ungz}) <(tail -n+2 {output.ungz}| sort -k2,2 -k3,3g) | bgzip -c > {output.gz} 
        tabix -s 2 -b 3 -e 3 {output.gz} 
        """


checkpoint polyAdapt_freqs:
    input:
        infile=os.path.join(config['uk_dir'], "{pheno}.flipped.byP.gz"),
        popfile="paneldir/{level}-euras.clusters.acf.gz",
        lbd=config['lbd']
    output:
        freqs="UKBiobank/data/gwasfreqs-{level}-{pheno}.tsv",
        outfile="UKBiobank/data/gwasfreqs-{level}-{pheno}.tsv.gz",
        candi="UKBiobank/data/gwasfreqs_candidates-{level}-{pheno}.tsv",
        neut="UKBiobank/data/gwasfreqs_neutral-{level}-{pheno}.tsv"
    shell:
        """
        python scripts/acf2ukbfreq_byP.py  -a {input.popfile} -g {input.infile} -o {output.freqs}
        cat <(head -1 {output.freqs}) <(tail -n+2 {output.freqs} | sort -k1,1 -k2,2g) | bgzip -c > {output.outfile}
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






