library(tidyr) ;
r <- dir() ;
s <- gsub("_S.*","",r) ;
a <- 0 ;
b <- 0 ;
b1 <- 0 ;
b2 <- 0 ;
for (i in 1:length(s)){ 
  a <- read.csv(r[i], header=T) 
  b <- s[i] 
  b1 <- append(b1,b)
  b2 <- append(b2,a)
}
b3 <- data.frame(sample=b1,data=as.character(b2))
b4 <- separate(b3,"data",into=c("reference","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq"),sep="\t") 
b5 <- b4[2:length(b4$reference),]
b5
write.table(b5,"coverages.tsv", sep="\t", row.names=F)
q("no") ;
