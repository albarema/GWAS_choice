import sys,os,pysam,gzip,random,numpy
from optparse import OptionParser
import subprocess

parser = OptionParser("$prog [options]")
parser.add_option("-i", "--infile", dest="infile", help="Input file (def None)", default=None, type="string")
parser.add_option("-o", "--outfile", dest="outfile", help="Output file (def None)", default=None, type="string")
(options,args) = parser.parse_args()


infile = gzip.open(options.infile,"r")


KGfiles = {}
for i in range(1,23):
    KGfiles[i] = pysam.Tabixfile("/willerslev/users-shared/science-snm-willerslev-rgf352/fernando/1000genomes/ALL.chr"+str(i)+".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",mode='r')

outfile = open(options.outfile,"w")

header = infile.readline().strip()
#header = header.replace("A1","REF")
#header = header.replace("A2","ALT")

print >>outfile, header
# a1 should be the effect size, most cases it is polarized by minor allele (a1 = minor). Now a1 = fields-4 for ukbv2. The other a1 = 3

for line in infile:

    fields = line.strip().split("\t")
    snpid = fields[0]
    chr = fields[1]
    pos = int(fields[2])
    a1 = fields[4]
    a2 = fields[3]
    Z = float(fields[5])
    SE = float(fields[6])
    P = float(fields[7])
    #print a1,a2
    elem = ""
    for extract in KGfiles[int(chr)].fetch(chr,pos-1,pos):
        elem = extract
    if elem == "": continue
    elemfields = elem.strip().split("\t")
    ref = elemfields[3]
    alt = elemfields[4]
    #print ref,alt

    # In GIANT data, effect allele is a1!
    if ref == a1 and alt == a2:
        Z = (-1)*Z
        print >>outfile, "\t".join([snpid,chr,str(pos),str(ref),str(alt),str(Z),str(SE),str(P)])
    elif ref == a2 and alt == a1:
        print >>outfile, "\t".join([snpid,chr,str(pos),str(ref),str(alt),str(Z),str(SE),str(P)])

        
    #print >>outfile, " ".join([Segnumber,SNPID,Chrom,Pos,PPA,Z,V])

