# This file performs the final editing of the Metabarpark_18S database.
# Sample counts are corrected by subtracting the mean counts in the blanks
# Then only MOTUs with ID> 0.80 are retained
# Then there is a low-frequencies filtering, which deletes counts of samples which are < 3% of the cumulative freq of each MOTU.
# Then, MOTUs with < 3 total counts are deleted
# Finally, incorrectly identified MOTUs with less than 200 counts are also deleted
# And incorrectly identified MOTUs with more than 200 counts are saved in another file for further manual checking

#load database
setwd("/media/owen/73E2-0BC8/metabarpark_18S") #Write your own working directory here
db <- read.csv("MBPARK18S.crop.count.good.blanks.tag.csv",as.is =T,head=T,sep="\t") #Select the correct separator. "\t" means tabulation.

#Add a new column with the initial sum of sample counts (sample columns are columns 30-108 in this file)
db$initial.counts <- rowSums(db[30:108])

#Calculate the means of the three blanks plus the negative (and round them to the next integer using ceiling)
blanks <- db[,26:29] #26 to 29 are the columns for the blanks and the "negative"
medias<-ceiling(rowMeans(blanks))

#Substract the mean of the blanks from all samples
db.corrected <- db[,30:108]-medias #30-108 are the columns for the samples.
db.corrected[db.corrected<0] <- 0 # All the negative values become 0
db.new <- data.frame(db[1:25],db.corrected,db[110:113]) #Remake the database, exchanging the sample columns for the corrected ones.

#Simplify the names of the samples
names(db.new)[26:104] <- substr(names(db.new)[26:104],8,nchar((db.new)[26:104]))

#Extract only MOTUs with best_ID>0.8
db.new <- db.new[db.new$best_identity.db_18s_r117>0.8,]

#Delete MOTUs with 0 counts
suma.MOTUs <- rowSums(db.new[26:104])
db.new <- db.new[suma.MOTUs>0,]

#Extract counts database
counts.db <- db.new[,26:104]
suma.samples <- colSums(counts.db,na.rm=T)

#Calculate frequency matrix
freq.matrix <- counts.db #For dimensioning
for (columna in 1:ncol(freq.matrix)) freq.matrix[,columna] <- counts.db[,columna]/suma.samples[columna]
check <- colSums(freq.matrix,na.rm=T)

#Calculate renormalized freq matrix
suma.taxa <- rowSums(freq.matrix)
freq.matrix.normalized <- freq.matrix #For dimensioning
for (fila in 1:nrow(freq.matrix)) freq.matrix.normalized[fila,] <- freq.matrix[fila,]/suma.taxa[fila]
check2 <- rowSums(freq.matrix.normalized)

#Delete all counts below 97% of cumulative freq
counts.trimmed <- counts.db
cutoff=0.97
for (fila in 1:nrow(freq.matrix)){
  #This orders the counts for each MOTU across all samples and calculates the cumulative value
  numeros <- as.vector(freq.matrix.normalized[fila,order(freq.matrix.normalized[fila,],decreasing=T)],"numeric")
  acumulado <- cumsum(numeros)
  muestra.limit <- sum(acumulado<cutoff) #Calculates the position of the cutoff value in the ordered vector   
  if (muestra.limit==0) muestra.limit <- 1 #If the position is the first, then there will be only 1 (and not 0).
  threshold <- numeros[muestra.limit] #This is the threshold frequency that will be different from 0
  #Now, all counts with a frequency below the threshold frequency of the MOTU, will be equated to 0
  for (columna in 1:ncol(freq.matrix)){
    if (freq.matrix.normalized[fila,columna]< threshold) counts.trimmed[fila,columna] <- 0
  }
  print(paste("fila:",fila))  #This is for monitoring the progress
}

#Reconstruct the database
new.db <- data.frame(db.new[,1:21],counts.trimmed,initial_counts=db.new[,108],
                     final_counts=rowSums(counts.trimmed),final_counts_Cabrera=rowSums(counts.trimmed[,1:36]),
                     final_counts_Cies=rowSums(counts.trimmed[,37:79]),sequence=db.new[,107])                     

#Delete all MOTUs with less than 3 counts
new.db <- new.db[new.db$final_counts>2,]

#Delete undesirable counts from sample CIE16A
new.db[(new.db$final_counts_Cies/new.db$final_counts_Cabrera > 0) & (new.db$final_counts_Cies/new.db$final_counts_Cabrera < 0.01)& (new.db$final_counts_Cies==new.db$Cie16A),c(68,104)] <- 0

#Correct the final counts
new.db$final_counts=rowSums(new.db[22:100])

#Add a column with trimmed counts
new.db$trimmed <- new.db$initial_counts-new.db$final_counts

#List samples to look by hand (unidentified with counts> 200)
look.by.hand <- new.db[new.db$scientific_name %in% c("Eukaryota","Eumetazoa", "Metazoa","Protostomia","Bilateria","Opisthokonta", "Ecdysozoa", "Lophotrochozoa", "Fungi", "Stramenopiles") & new.db$final_counts > 200,]
write.table(look.by.hand,"Metabarpark_18S_to_look_by_hand.csv",sep=";",row.names=F)

#Delete unidentified sampled with low counts
borrar <- new.db[new.db$scientific_name %in% c("Eukaryota","Eumetazoa", "Metazoa","Protostomia","Bilateria","Opisthokonta", "Ecdysozoa", "Lophotrochozoa", "Fungi", "Stramenopiles") & new.db$final_counts <= 200,]
write.table(borrar,"Metabarpark_18S_auto_deleted.csv",sep=";",row.names=F)
new.db <- new.db[!(new.db$scientific_name %in% c("Eukaryota","Eumetazoa", "Metazoa","Protostomia","Bilateria","Opisthokonta", "Ecdysozoa", "Lophotrochozoa", "Fungi", "Stramenopiles") & new.db$final_counts <= 200),]

#Save the final database in csv format
write.table(new.db,"Metabarpark_18S_new.db.csv",sep=";",row.names=F)



