#!/usr/bin/env Rscript
library("optparse")
# terminal line
# Rscript CalcQX.R -w gwasfreqs_candidates_pheno.tsv -e gwasfreqs_neutral_pheno.tsv -o qxreport.txt -n 1000 -j 1000 -s genscores.txt 

option_list = list(
    make_option(c("-w", "--gwasfile"), type="character", default=NULL, help="GWAS input file name"),
    make_option(c("-e", "--neutfile"), type="character", default=NULL, help="Neutral input file name"),
    make_option(c("-o", "--outfile"), type="character", default="output.txt", help="Output file name"),
    make_option(c("-p", "--pops"), type="character", default="ALL", help="Populations to be tested"),
    make_option(c("-n", "--numrep"), type="numeric", default=NULL, help="Number of sign-randomized replicates"),
    make_option(c("-j", "--emprep"), type="numeric", default=NULL, help="Number of frequency-matched replicates"),
    make_option(c("-s", "--scorefile"), type="character", default=NULL, help="Score file name"),
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
if (is.null(opt$outfile)){
    print <- help(opt_parser)
    stop("Output file name must be supplied.n", call.=FALSE)
}

if (is.null(opt$scorefile)){
    print <- help(opt_parser)
    stop("Score file name must be supplied.n", call.=FALSE)
}


gwasfile <- opt$gwasfile
neutfile <- opt$neutfile
outfile <- opt$outfile
scorefile <- opt$scorefile
pops <- opt$pops
#pops <- strsplit(opt$pops,",")[[1]]
pseudorep <- opt$numrep
emprep <- opt$emprep

# Multiple-testing
popmatch <- opt$popmatch
leeway <- opt$leeway
minorfreq <- opt$minorfreq


#pops <- strsplit(readLines(pops), ",")[[1]]
#pops <- c("ESN","MSL", "CEU", "TSI", "CDX", "JPT", "PEL")

ObtainFreqs <- function(countdat){
    dercounts <- apply(countdat,c(1,2),function(x){splitted <- strsplit(x,",")[[1]]; return( as.numeric(splitted[2]) )})
    totalcounts <- apply(countdat,c(1,2),function(x){splitted <- strsplit(x,",")[[1]]; return( as.numeric(splitted[2])+as.numeric(splitted[1]) )})
    dersum <- apply(dercounts,1,function(x){sum(x)})
    totalsum <- apply(totalcounts,1,function(x){sum(x)})
    totalfreq <- dersum/totalsum
    freqs <- apply(countdat,c(1,2),function(x){splitted <- strsplit(x,",")[[1]];
        return( as.numeric(splitted[2]) / (as.numeric(splitted[2])+as.numeric(splitted[1])) )})
    return(list(as.matrix(freqs), totalfreq, as.matrix(dercounts),as.matrix(totalcounts)))
}

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

ComputeFmat <- function(neut_leaves_freqs, neut_total_freqs){
    leaves <- colnames(neut_leaves_freqs)

    #checksegneut <- which( apply(neut_leaves_freqs,1,sum)/length(leaves) < 0.95  & apply(neut_leaves_freqs,1,sum)/length(leaves) > 0.05 )
    #neut_leaves_freqs <- neut_leaves_freqs[checksegneut,]
    checksegneut <- which(neut_total_freqs < 0.95  & neut_total_freqs > 0.05 )
    neut_leaves_freqs <- neut_leaves_freqs[checksegneut,]
    
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


ChiSquared <- function(leaves_freqs,total_freqs,effects,F_mat, dercounts, totalcounts, randomize=FALSE){
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
    
    # Compute error PRS
    alpha <- dercounts + 1
	beta <- 1 + totalcounts - dercounts
	num = alpha * beta
	denom = (alpha+beta)^2 * (alpha+beta+1)
	varpl = num/denom

	varpl <- varpl[c(checkseg),]
	varprs = effects^2 * varpl

	varprs <- apply(varprs, 2, function(x) sum(x))
	se = sqrt(varprs)


	# meangen already subracted mean 
	normcilow <- 2*(meangen-1.96*se)/sqrt(4*varmean)
	normcihi <- 2*(meangen + 1.96*se)/sqrt(4*varmean)

    # Compute Q_X statistic
    numerator <- t(meangen) %*% solve(Fmat) %*% meangen
    denominator <- varmean
    Qteststat <- numerator / denominator
    Pval <- 1 - pchisq(Qteststat,qr(Fmat)$rank)
    allstats <- c(Qteststat,Pval)

    return(list(allstats, meangenvec, normcihi, normcilow))
}


# Load GWAS data  
print("Loading data...")
data <- LoadCounts(gwasfile, pops)
leaves_counts <- as.data.frame(data[,seq(5,dim(data)[2])])
raw_leaves_freqs <- ObtainFreqs(leaves_counts)
leaves_freqs <- raw_leaves_freqs[[1]]
total_freqs <- raw_leaves_freqs[[2]]
dercounts <- raw_leaves_freqs[[3]]
totalcounts <- raw_leaves_freqs[[4]]
effects <- as.numeric(data[,4])


if (sum(apply(leaves_freqs, 1, function(x) sum(is.na(x))))!=0){
leaves_freqs <- na.omit(leaves_freqs)
total_freqs <- na.omit(total_freqs)
print(paste('Ommiting', length(leaves_counts[,1])-length(leaves_freqs[,1]),
            'out of', length(leaves_counts[,1]), 'character SNPs', sep = ' '))
}

# Load Neutral data
neutdata <- LoadCounts(neutfile, pops)
neut_leaves_counts <- as.data.frame(neutdata[,seq(5,dim(neutdata)[2])])
raw_neut_leaves_freqs <- ObtainFreqs(neut_leaves_counts)
neut_leaves_freqs <- raw_neut_leaves_freqs[[1]]
neut_total_freqs <- raw_neut_leaves_freqs[[2]]

# Calculate covariance matrix
print("Computing covariance matrix...")
Fmat <- ComputeFmat(neut_leaves_freqs, neut_total_freqs)


# Calculate chi-squared statistics
print("Computing Q_X statistic...")
totaltest <- ChiSquared(leaves_freqs,total_freqs,effects,Fmat, dercounts, totalcounts, randomize=FALSE)
totalstat <- totaltest[[1]]
meangenvec <- totaltest[[2]]
meangenvechi <- totaltest[[3]]
meangenveclow <- totaltest[[4]]


qtab <- cbind(round(totalstat[1],3),totalstat[2])
colnames(qtab) <- c("Q_X","Pval")

#### Only calculate pval and genetic scores

write("Genetic scores",file=outfile,append=FALSE)
write(paste(names(meangenvec),collapse="\t"),file=outfile,append=TRUE)
write(paste(meangenvec,collapse="\t"),file=outfile,append=TRUE)
write("",file=outfile,append=TRUE)

printgenvec <- cbind(names(meangenvec),meangenvec, meangenvechi, meangenveclow)
rownames(printgenvec) <- c()
colnames(printgenvec) <- c("#POP","SCORE", "upperlim", "lowerlim")
write.table(printgenvec,file=scorefile,sep="\t",row.names=FALSE,col.names=TRUE, quote=FALSE,append=FALSE)

write(paste("Q_X:",qtab[1],sep="\t"),file=outfile,append=TRUE)
write("",file=outfile,append=TRUE)

write(paste("Chi-squared distribution P-value:",qtab[2],sep="\t"),file=outfile,append=TRUE)
write("",file=outfile,append=TRUE)

# Add number of traits-ass SNPs 
write(paste("Number of trait-associated SNPs:", length(effects), sep="\t"), file=outfile, append=TRUE)
write("",file=outfile,append=TRUE)

# Calculate sign-randomized chi-squared statistics
if( !is.null(pseudorep) ){
    print(paste("Computing sign-randomized P-values, using ",pseudorep," pseudo-replicates...",sep=""))
    totaltest <- ChiSquared(leaves_freqs,total_freqs,effects,Fmat, dercounts, totalcounts, randomize=TRUE)
    totalstat <- totaltest[[1]]
    randqtab <- round(totalstat[1],3)
    for(i in seq(2,pseudorep)){
        totaltest <- ChiSquared(leaves_freqs,total_freqs,effects,Fmat, dercounts, totalcounts, randomize=TRUE)
        totalstat <- totaltest[[1]]
        randqtab <- c( randqtab, round(totalstat[1],3) )
    }
    # Edit adding 1
    randpval <- (1+sum( as.numeric(randqtab)  > as.numeric(qtab[1]) )) / (1+length(randqtab))
    
    write(paste("Sign-randomized P-value:",randpval,sep="\t"),file=outfile,append=TRUE)
    write("",file=outfile,append=TRUE)
    
}

