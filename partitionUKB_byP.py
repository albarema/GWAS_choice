import pysam
import numpy as np
from optparse import OptionParser
import subprocess

parser = OptionParser("$prog [options]")
parser.add_option("-i", "--infile", dest="infile", help="Input GWAS+freq file", default=None, type="string")
parser.add_option("-b", "--bedfile", dest="bedfile", help="Bed file with LD partitions (def None)", default=None, type="string")
parser.add_option("-p", "--maxp", dest="maxp", help="Maximum p-value allowed", default=1, type="float")
parser.add_option("-o", "--outfile", dest="outfile", help="Output file", default=None, type="string")
(options,args) = parser.parse_args()


bedfile = open(options.bedfile,"r")
infile = pysam.Tabixfile(options.infile, mode='r')
outfile = open(options.outfile,"w")
logmaxp = np.log10(options.maxp)

readheader = infile.header
for x in readheader:
    header = x

header = header.split("\t")
header = list(filter(lambda a: a != "REF" and a != "ALT" and a != "ANC" and a != "DER" and a!= "SE" and a != "PVAL", header))
header = "\t".join(header)
outfile.write(''.join(header) + '\n')

for line in bedfile:
    regfields = line.strip("\n").split("\t")
    regchr = regfields[0].strip()
    regstart = int(regfields[1])
    regend = int(regfields[2])

    CurrLogPval = 0
    CurrPval = None
    CurrBest = None
    try:
        for elem in infile.fetch(regchr,regstart-1,regend):
            fields = elem.strip().split("\t")
            Pval = float(fields[9])
            if Pval == 0.0:
                Pval = 1e-20
            LogPval = np.log10(Pval)
            if LogPval < CurrLogPval:
                CurrBest = fields[0:3]+[fields[7]]+fields[10:]
                CurrLogPval = LogPval
                CurrPval = Pval
    except:
        continue

    if CurrLogPval < logmaxp and CurrBest != None:
        outfile.write("\t".join(CurrBest) + '\n')


