#!/usr/bin/python
#!/usr/bin/python

import sys,os,pysam,itertools,random,numpy,gzip, csv
from optparse import OptionParser

parser = OptionParser("$prog [options]")
parser.add_option("-a", "--acffile", dest="acffile", help="Input ACF file", default=None, type="string")
parser.add_option("-o", "--outfile", dest="outfile", help="Output file (def None)", default=None, type="string")
parser.add_option("-g", "--gwasfile", dest="gwasfile", help="GWAS file", default=None, type="string")
parser.add_option("-l", "--ldblockfile", dest="ldblockfile", help="LD block file (def None)", default=None, type="string")
(options,args) = parser.parse_args()


# Read vcf and ancestral state files
acffile = pysam.Tabixfile(options.acffile,mode='r')

# Open output file
outfile = open(options.outfile,"w")

# Record population panels
for line in acffile.header:
    if "coor" in line:
        popordered = line.split("\t")[5:]

# Print header
header = "#CHROM\tPOS\tSNPID\tREF\tALT\tANC\tDER\tDEREFFECT\tSE\tPVAL\t"
header += "\t".join(popordered)
outfile.write(''.join(header) + '\n')


# Read files
with gzip.open(options.gwasfile, mode="rt") as tsvfile:
    reader = csv.DictReader(tsvfile, dialect='excel-tab')
    for row in reader:
        info = row['variant']
        infofields = info.strip("\n").split(":")
        chrom = infofields[0]
        pos = infofields[1]
        gref = infofields[2]
        galt = infofields[3]
        effect = row['beta']
        se = row['se']
        pval = row['pval']
        snpid = str(chrom) + ":" + str(pos)
        prevpos = int(pos) - 1
        elem = "NA"
        try:
            acfline = acffile.fetch(str(chrom),int(prevpos),int(pos))
        except:
            continue
        for subline in acfline:
            elem = subline
            if elem == "NA":
                continue
            fields = elem.strip("\n").split("\t")
            refalt = fields[2].split(",")
            ref = refalt[0]
            alt = refalt[1]
            ancinfo = fields[4].split(":")[0]

            if ancinfo == "1,0":
                anc = ref
                der = alt
            elif ancinfo == "0,1":
                anc = alt
                der = ref
            else:
                continue

            if alt == galt:
                if anc != ref:
                    effect = str((-1)*float(effect))
                    popfreqs = [",".join(list(reversed(x.split(":")[0].split(",")))) for x in fields[5:]]
                else:
                    popfreqs = [x.split(":")[0] for x in fields[5:]]
                finalvec = [chrom,pos,snpid,gref,galt,anc,der,effect,se,pval]+popfreqs
            else:
                continue

            outfile.write("\t".join(finalvec) + '\n')

outfile.close()
