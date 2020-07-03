import pysam
from optparse import OptionParser

parser = OptionParser("$prog [options]")
parser.add_option("-i", "--infile", dest="infile", help="Input GWAS+freq file", default=None, type="string")
parser.add_option("-p", "--minp", dest="minp", help="Minimum p-value allowed", default=1, type="float")
parser.add_option("-o", "--outfile", dest="outfile", help="Output file", default=None, type="string")
parser.add_option("-s", "--sep", dest="sep", help="Number of valid SNPs separating each printed SNP", default=None, type="int")
(options,args) = parser.parse_args()


infile = pysam.Tabixfile(options.infile, mode='r')
outfile = open(options.outfile,"w")

readheader = infile.header
for x in readheader:
    header = x

header = header.split("\t")
header = list(filter(lambda a: a != "REF" and a != "ALT" and a != "ANC" and a != "DER" and a!= "SE" and a != "PVAL", header))
header = "\t".join(header)

outfile.write(''.join(header) + '\n')

i=0
for line in infile.fetch():
    fields = line.strip().split("\t")
    Pval = float(fields[9])
    if Pval > options.minp:
        if i == options.sep:
            i = 0
            finalvec = fields[0:3]+[fields[7]]+fields[10:]
            outfile.write("\t".join(finalvec)+ '\n')
        i += 1


