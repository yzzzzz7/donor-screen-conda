args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
p <- args[1]
library(RZooRoH)
stopifnot(as.character(packageVersion("RZooRoH"))=="0.4.3")
rates <- c(5,25,125,625,3125,15625)
model <- zoomodel(K=6,predefined=TRUE,krates=rates,err=0.001,seqerr=0.001)
ids <- readLines(paste0(p,".samples"))
freq <- scan(paste0(p,".freq"),quiet=TRUE)
data <- zoodata(paste0(p,".haps"),zformat="haps",samplefile=paste0(p,".samples"),
                allelefreq=freq,min_maf=0,map_unit="bp",RecRate=1)
stopifnot(identical(as.character(data@sample_ids),ids),data@nsnps==length(freq))
pairs <- as.matrix(read.delim(paste0(p,".pairs")))
answer <- numeric(nrow(pairs))
for (start in seq.int(1L,nrow(pairs),by=8L)) {
    i <- seq.int(start,min(start+7L,nrow(pairs)))
    kin <- pairs[i,,drop=FALSE]
    fit <- zookin(model,data,kinpairs=kin,nT=128,fb=TRUE,vit=FALSE,localhbd=FALSE)
    stopifnot(is(fit,"kres"),all(dim(fit@realized)==c(length(i),6)),
              identical(as.character(fit@sampleids),paste(ids[kin[,1]],ids[kin[,2]],sep="_")),
              length(fit@optimerr)==length(i),all(fit@optimerr==0),
              all(is.finite(fit@realized)),all(fit@realized>=0),
              all(rowSums(fit@realized)<=1+1e-8))
    saveRDS(fit,sprintf("%s.batch_%04d.rds",p,start))
    answer[i] <- rowSums(fit@realized[,1:2,drop=FALSE])
    cat("Completed pairs",max(i),"/",nrow(pairs),"\n")
    flush.console()
}
writeLines(sprintf("%.17g",answer),paste0(p,".values"))

