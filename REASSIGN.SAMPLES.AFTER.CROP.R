# This script will join information from three sources:
# - A tabulated file from obitools including sequence counts for individual samples.
# - A list of MOTU clusters joined by CROP (with the same ids than the previous file)
# - A file including fasta sequences for cluster-heads from a cluster.fasta.file and taxonomic info 
#   for these sequences obtained by ecotag (with the same ids than the previous files)
# All this will be joined in a csv output file that may be opened with Calc, Excel or R for further calculations
# The counts of every sequence in the first file will be added to the cluster-head where the sequence belongs,
# keeping the total counts and the counts for each sample.

setwd("~/metabarpark_18S")
filetab <- "MBPARK18S.tab.txt" #Must be a tabulated file (from obitab) including individual sample counts
fileclust <- "MBPARK18S_1.filtered.noSingleton.len.fasta.cluster.list" #Must be a cluster-list outfile from CROP
filetag <- "MBPARK18S.crop.tag.tab.txt" #Must be a taxonomy-annotated file from applying ecotag to a cluster.fasta file and then tabulated with obitab
outfile <- "MBPARK18S.crop.count.tag.csv"

db <- read.table(filetab,sep="\t",head=T)
db$id <- as.character(db$id)
db <- db[-c(2,4,5,7,9:12,14,107,108,110:125,127:129)]
numseqs <- nrow(db)
clusters <- read.table(fileclust,sep="\t",head=F,stringsAsFactors=F)
numclust <- nrow(clusters)
ecotag <- read.table(filetag,sep="\t",head=T,stringsAsFactors=F)
ecotag <- ecotag[-c(2,5)]
db.total <- merge(ecotag,db,by="id") #This will keep just the heads
db.total <- data.frame(db.total[-c(ncol(db.total))],cluster_count = 1,db.total$sequence)
clusters.ord <- data.frame(id=clusters$V1[order(clusters$V1)],included=clusters$V2[order(clusters$V1)])
rm(clusters)
for (fila in 1:numclust){
  head <- clusters.ord[fila,1]
  tails <- strsplit(as.character(clusters.ord[fila,2]),",")
  for (j in 1:length(tails[[1]])){
    tail <- tails[[1]][j]
    if (tail!=head) { 
      db.total[fila,c(22,26:116)] <- db.total[fila,c(22,26:116)] + db[db$id==tail,c(3,7:97)] 
      db.total[fila,119] <- db.total[fila,119]+1
    } 
  }
  print(paste("Cluster", fila, "/",numclust, "ready"))
}
write.table(db.total,outfile,sep=";",quote=F,row.names=F)
print(paste("File ", outfile, " written"))
