# use this code to install sequenza in r:
# setRepositories(graphics = FALSE, ind = 1:6)
# install.packages("sequenza")
library("sequenza")
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
data.file <- args[2]
seqzdata <- sequenza.extract(data.file)
CP.example <- sequenza.fit(seqzdata)
sequenza.results(sequenza.extract = seqzdata, cp.table = CP.example, sample.id = args[3], out.dir = args[3])
