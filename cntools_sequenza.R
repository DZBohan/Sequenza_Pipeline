library(CNTools)

# set the arguments
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
gene.info.f <- args[2]
sqza_file <- args[3]
id <- args[4]

# read the geninfo file
geneinfo.hg38 = read.table(gene.info.f,header=T,sep="\t")

# read the results segment files of sequenza
sqza_out <- read.table(sqza_file,header=T)

# standardize the table of sequenza output
sqza_out <- sqza_out[,c(-4:-6,-9:-13)]
sqza_out$ID<-c(id)
sqza_out$seg.mean <- log(sqza_out$depth.ratio, 2)
sqza_out <- sqza_out[,c(-4)]
names(sqza_out)[1] <- "chrom"
names(sqza_out)[2] <- "loc.start"
names(sqza_out)[3] <- "loc.end"
names(sqza_out)[4] <- "num.mark"
sqza_out <- sqza_out[,c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")]# the below line is require when there are "chr" for the chrom colums
sqza_out$chrom <- substring(sqza_out$chrom, 4)
sqza_out <- na.omit(sqza_out)

# run cntools on the segments results of sequenza
cnseq_sqza <- CNSeg(sqza_out)
rd_sqza <- getRS(cnseq_sqza, by="gene", imput=FALSE, XY=FALSE, geneMap=geneinfo.hg38, what = "mean")
rs_sqza <- rs(rd_sqza)

# remove 0 value in the table 
rs_sqza[rs_sqza==0] <- NA
rs_sqza <- na.omit(rs_sqza)

# output the txt file
write.table (rs_sqza, file = paste(id,"cntools","sequenza","txt",sep = "."), quote =FALSE, sep ="\t", row.names=FALSE)
