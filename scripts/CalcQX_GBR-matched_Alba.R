#!/usr/bin/env Rscript
suppressWarnings(suppressMessages(library("optparse")))
suppressWarnings(suppressMessages(library(tidyverse)))
# terminal line
# Rscript scripts/CalcQX_GBR-matched_Alba.R -w {input.candi} -e {input.neut} -a {input.gbr} -n 1000 -m {output.qxfm} -j 1000
# input_gbr is a file in allele count format generated with glactools only targeting GBR population. Useful when we are working with a dataset that does not include 1000GP but we still want to use it

option_list = list(
    make_option(c("-w", "--gwasfile"), type="character", default=NULL, help="GWAS input file name"),
    make_option(c("-e", "--neutfile"), type="character", default=NULL, help="Neutral input file name"),
    make_option(c("-a", "--acffile"), type="character", default=NULL, help="1000g-pop count input file name"),
    make_option(c("-m", "--outfile2"), type="character", default="output2.txt", help="Output file name fmatch"),
    make_option(c("-p", "--pops"), type="character", default="ALL", help="Populations to be tested"),
    make_option(c("-n", "--numrep"), type="numeric", default=NULL, help="Number of sign-randomized replicates"),
    make_option(c("-j", "--emprep"), type="numeric", default=NULL, help="Number of frequency-matched replicates"),
    make_option(c("-l", "--leeway"), type="numeric", default=0.01, help="Leeway for empirical frequency matching"),
    make_option(c("-f", "--minorfreq"), type="numeric", default=0.05, help="Minor allele frequency cutoff for covariance matrix")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$gwasfile)){
    print <- help(opt_parser)
    stop("GWAS file name must be supplied.n", call.=FALSE)
}
if (is.null(opt$neutfile)){
    print <- help(opt_parser)
    stop("Neutral file name must be supplied.n", call.=FALSE)
}

if (is.null(opt$acffile)){
    print <- help(opt_parser)
    stop("Acf file name must be supplied.n", call.=FALSE)
}

if (is.null(opt$outfile2)){
    print <- help(opt_parser)
    stop("Output file name2 must be supplied.n", call.=FALSE)
}


gwasfile <- opt$gwasfile
neutfile <- opt$neutfile
acffile <- opt$acffile
outfile2 <- opt$outfile2
pseudorep <- opt$numrep
emprep <- opt$emprep

# Multiple-testing
popmatch <- opt$popmatch
leeway <- opt$leeway
minorfreq <- opt$minorfreq

# Uncomment this if not using ALL populations (otherwise pops are comma separated in separate file)
# pops <- strsplit(readLines(pops), ",")[[1]]
pops <- "ALL"

# Trim out white space
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

LoadCounts <- function(filename,pops){
    table <- as.matrix(read.table(filename,header=TRUE,sep="\t",strip.white=TRUE,comment.char="!"))
    if(paste(pops,collapse=",") == "ALL"){ tokeep <- seq(5,length(colnames(table)))
    } else{ tokeep <- which(colnames(table) %in% pops)}
    tokeep <- c(1,2,3,4,tokeep)
    table <- table[,tokeep]
    table[,1] <- trim(table[,1])
    table[,2] <- trim(table[,2])
    return(table)
}

ObtainFreqs <- function(countdat){
    dercounts <- apply(countdat,c(1,2),function(x){splitted <- strsplit(x,",")[[1]]; return( as.numeric(splitted[2]) )})
    totalcounts <- apply(countdat[,-1],c(1,2),function(x){splitted <- strsplit(x,",")[[1]]; return( as.numeric(splitted[2])+as.numeric(splitted[1]) )})
    dersum <- apply(dercounts,1,function(x){sum(x)})
    totalsum <- apply(totalcounts,1,function(x){sum(x)})
    totalfreq <- dersum/totalsum
    freqs <- apply(countdat,c(1,2),function(x){splitted <- strsplit(x,",")[[1]];
        return( as.numeric(splitted[2]) / (as.numeric(splitted[2])+as.numeric(splitted[1])) )})
    return(list(as.matrix(freqs), totalfreq))
}

ComputeFmat <- function(neut_leaves_freqs, neut_total_freqs){
    leaves <- colnames(neut_leaves_freqs)
    checksegneut <- which(neut_total_freqs < 0.95  & neut_total_freqs > 0.05 )
    neut_leaves_freqs <- neut_leaves_freqs[checksegneut,]
    d <- apply(neut_leaves_freqs, 1, sum)
    tokeep <- which(d != 0)
    neut_leaves_freqs <- neut_leaves_freqs[tokeep,]
    
    neut_leaves_freqs_means <- apply(neut_leaves_freqs, 1, mean)
    mean_hetero <- neut_leaves_freqs_means*(1-neut_leaves_freqs_means)
    numSNPs <- length(neut_leaves_freqs_means)

    Fmat <- sapply(seq(1,dim(neut_leaves_freqs)[2]),function(x){
        sapply(seq(1,dim(neut_leaves_freqs)[2]),function(y){
            cov(neut_leaves_freqs[,x]/sqrt(mean_hetero), neut_leaves_freqs[,y]/sqrt(mean_hetero))
        })
    })

    colnames(Fmat) <- colnames(neut_leaves_freqs)
    rownames(Fmat) <- colnames(neut_leaves_freqs)
    return(Fmat)
}


ChiSquared <- function(leaves_freqs,total_freqs,effects,F_mat,randomize=FALSE){
    leaves <- colnames(leaves_freqs)
    checkseg <- which(total_freqs < 0.95  & total_freqs > 0.05 )
    leaves_freqs <- leaves_freqs[c(checkseg),]
    effects <- effects[c(checkseg)]
    leaves_freqs <- leaves_freqs[! is.infinite(effects),]
    effects <- effects[! is.infinite(effects)]

    # Randomize effects if necessary
    if(randomize == TRUE){effects <- effects * sample(c(-1,1),length(effects),replace=TRUE)}
    
    # Compute mean genetic values
    meangen <- apply(leaves_freqs * effects, 2, function(x){sum(x)})
    # Scale by average genetic value
    meangen <- (meangen - mean(meangen))

    # Compute the estimated ancestral genetic variance over all populations
    meanfreqs <- apply(leaves_freqs,1,mean)

    varmean <- sum(sapply(seq(1,length(meanfreqs)),function(i){
        score = meanfreqs[i]*(1-meanfreqs[i])*effects[i]^2
        return(score)
    }))

    meangenvec <- 2*meangen / sqrt( 4*varmean)
    
    # Compute Q_X statistic
    numerator <- t(meangen) %*% solve(Fmat) %*% meangen
    denominator <- varmean
    Qteststat <- numerator / denominator
    Pval <- 1 - pchisq(Qteststat,qr(Fmat)$rank)
    if(Pval < 2.220446e-16){ Pval = 2.220446e-16}
    allstats <- c(Qteststat,Pval)

    return(list(allstats,meangenvec))
}

New_can <- function(acf, gwas){
	gwas <- gwas[,-c(3,4)]
	cand <- merge(acf, gwas, by=c(1,2))
}

New_neut <- function(acf, neut){
	neut <- neut[,-c(3,4)]
	neut <- merge(acf,neut, by=c(1,2))
}


# Load GWAS data  
print("Loading data...")
d <- LoadCounts(gwasfile, pops)
# Load Neutral data
neutd <- LoadCounts(neutfile, pops)
# Load ACF file and merge GBR with other pops
acff <- read_tsv(acffile)
data <- New_can(acff, d)
neutdata <- New_neut(acff, neutd)

# Get frequencies
#1. GWAS data
leaves_counts <- as.data.frame(data[,seq(3,dim(data)[2])])
raw_leaves_freqs <- ObtainFreqs(leaves_counts)
leaves_freqs <- raw_leaves_freqs[[1]]
total_freqs <- raw_leaves_freqs[[2]]
effects <- as.numeric(d[,4])


if (sum(apply(leaves_freqs, 1, function(x) sum(is.na(x))))!=0){
leaves_freqs <- na.omit(leaves_freqs)
total_freqs <- na.omit(total_freqs)
print(paste('Ommiting', length(leaves_counts[,1])-length(leaves_freqs[,1]),
            'out of', length(leaves_counts[,1]), 'character SNPs', sep = ' '))
}
#2. Neutral data
neut_leaves_counts <- as.data.frame(neutdata[,seq(3,dim(neutdata)[2])])
neutral_leaves_counts <- neut_leaves_counts %>% filter_all(all_vars( .!= "0,0"))
raw_neut_leaves_freqs <- ObtainFreqs(neutral_leaves_counts)
neut_leaves_freqs <- raw_neut_leaves_freqs[[1]]
neut_total_freqs <- raw_neut_leaves_freqs[[2]]

# Calculate covariance matrix
print("Computing covariance matrix...")
Fmat <- ComputeFmat(neut_leaves_freqs[,-1], neut_total_freqs)

# Calculate chi-squared statistics
print("Computing Q_X statistic...")
totaltest <- ChiSquared(leaves_freqs[,-1],total_freqs,effects,Fmat,randomize=FALSE)
totalstat <- totaltest[[1]]
meangenvec <- totaltest[[2]]


qtab <- cbind(round(totalstat[1],3),totalstat[2])
colnames(qtab) <- c("Q_X","Pval")


# Calculate frequency-matched chi-squared statistics
popmatch <- 'GBR'
for(i in popmatch){

    if( !is.null(i) & !is.null(emprep) ){
        print(paste("Computing P-values from frequency-matched empirical distribution, using ",emprep," ",i,"-matched replicates...",sep=""))
        allsamplestats <- sapply(seq(1,emprep),function(rep){

            sampleidx <- sapply(seq(1,dim(leaves_freqs)[1]),function(testsnp){
                sampleset <- which( neut_leaves_freqs[,i] > leaves_freqs[testsnp,i]-leeway & neut_leaves_freqs[,i] < leaves_freqs[testsnp,i]+leeway )
                snpidxsample <- sample( sampleset ,1)
                return(snpidxsample)
            })

            sample_leaves_freqs <-  neut_leaves_freqs[sampleidx,]
            totaltest <- ChiSquared(sample_leaves_freqs[,-1],total_freqs,effects,neut_leaves_freqs[,-1],randomize=FALSE)
            totalstat <- totaltest[[1]]
            qstat <- totalstat[1]
            return(qstat)
        })
    }
    # Edit with adding 1
    emppval <- (1+sum( as.numeric( allsamplestats )  > as.numeric(qtab[1]) )) / (1+length( allsamplestats ))

    write(paste("GBR-Frequency-matched P-value:", emppval,sep="\t"),file=outfile2,append=TRUE)

}

