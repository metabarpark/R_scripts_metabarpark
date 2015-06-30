## Script for simplifying an obitab output file by deleting everything but sequence ids and numbers of reads for every sample
## The script will look for columns beginning with "sample:" and will keep only the numbers under these columns
## The script will read two arguments from the command line: the input file name and the output file name.
## If the second argument is empty, it will just add ".csv" at the end of the name of the input file.
## By Owen S. Wangensteen - Project Metabarpark  2015

args<-commandArgs(T)

if (length(args)==0) {message("Please enter a tabulated input file you want to simplify.")} else
{  
filetab <- args[1] #Must be a tabulated file (from obitab) including individual sample counts under "sample:" columns
if (length(args)<2) outfile <-paste(filetab,".csv",sep="") else outfile  <- args[2] #If args[2] is empty, then  

entrada <- file(filetab,"r")
salida <- file(outfile,"w")

names <- as.vector(strsplit(readLines(entrada,1),"\t"))
newnames <- names[[1]][lapply(names,substr,start=1,stop=7)[[1]]=="sample:" | lapply(names,substr,start=1,stop=7)[[1]]=="id"]
writeLines(sub("sample:","",newnames[1:length(newnames)-1]),salida,sep=",")
writeLines(sub("sample:","",newnames[length(newnames)]),salida)
numberlines <- 0
repeat {
  linea <- as.vector(strsplit(readLines(entrada,1),"\t"))
  if (length(linea)==0) break
  newline <- linea[[1]][lapply(names,substr,start=1,stop=7)[[1]]=="sample:" | lapply(names,substr,start=1,stop=7)[[1]]=="id"]
  writeLines(newline[1:length(newline)-1],salida,sep=",")
  writeLines(newline[length(newline)],salida)
  numberlines <- numberlines+1
  if (numberlines %% 50 == 0) message(numberlines," lines written.","\r",appendLF = FALSE)
}
close(entrada)
close(salida)
message("Output file ",outfile," written with ",numberlines, " sequences")
}
