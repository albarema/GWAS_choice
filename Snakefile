#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:01:16 2020
@author: gsd818
"""


import pandas as pd
import csv

configfile: "config.yaml"

## --------------------------------------------------------------------------------
##### Modules #####
include: "rules/gwasApproaches.smk"
include: "rules/polyadapt.smk"
include: "vcf2acf.smk"

## --------------------------------------------------------------------------------
##### Wildcards #####
wildcard_constraints:
    pheno="[^-]+"
        
## --------------------------------------------------------------------------------
##### Target rules #####
def run_all_rules(_):
    inputs = []


    # this will run one of the phenotypes
    #"UKBiobank/data/gwasfreqs-pops-102_irnt.tsv.gz",

    this will run all the phenotypes
    with open("phenoname_qx.txt", "w") as fout:
        for pheno in pd.read_table('phenoname.txt')['phenoname'].tolist():
            tsv = checkpoints.polyAdapt_freqs.get(pheno=pheno, level='pops').output.candi

            with open(tsv) as fin:
                if len(fin.readlines()) > 4:
                    inputs.append("UKBiobank/selection_UKBV2/Genscores-pops-{pheno}.txt".format(pheno=pheno))
                    inputs.append("UKBiobank/selection_UKBV2/QX_fm_report-pops-{pheno}.txt".format(pheno=pheno))
                    
                    
## targets
rule allrules:
    input:
        run_all_rules
