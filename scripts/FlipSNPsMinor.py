#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 12:59:19 2020

@author: gsd818
"""

import csv, gzip, pysam
from optparse import OptionParser

parser = OptionParser("$prog [options]")
parser.add_option("-i", "--infile", dest="infile", help="Input file (def None)", default=None, type="string")
parser.add_option("-o", "--outfile", dest="outfile", help="Output file (def None)", default=None, type="string")
(options,args) = parser.parse_args()



KGfiles = {}
for i in range(1,23):
    KGfiles[i] = pysam.Tabixfile("/willerslev/datasets/1000_Genomes_phase3_v5a/individual_chromosomes/chr"+str(i)+".1kg.phase3.v5a.vcf.gz",mode='r')

outfile = open(options.outfile,"w")

header = "#SNPID\tCHR\tPOS\tREF\tALT\tEFFECT\tSE\tPVAL"
outfile.write(''.join(header) + '\n')

with gzip.open(options.infile,"r") as infile:
    #infile.readline()
    for line in infile:
        fields = line.decode('utf8').strip("\n").split()
        snpid = fields[2]
        chr = fields[0]
        pos = int(fields[1])
        a1 = fields[5]
        try:
            Z = float(fields[8])
            SE = float(fields[9])
            P = float(fields[11])
        except:
            continue
        #print a1,a2
        elem = ""
        for extract in KGfiles[int(chr)].fetch(chr,pos-1,pos):
            elem = extract
            if elem == "": continue
            elemfields = elem.strip().split("\t")
            ref = elemfields[3]
            alt = elemfields[4]
            # In UKB data, effect allele is a1!
        if ref == a1:
            Z = (-1)*Z
            outfile.write("\t".join([snpid,chr,str(pos),str(ref),str(alt),str(Z),str(SE),str(P)])+ '\n')
        elif alt == a1:
            outfile.write("\t".join([snpid,chr,str(pos),str(ref),str(alt),str(Z),str(SE),str(P)])+ '\n')


outfile.close()
