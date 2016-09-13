#Cat a bunch of files to their matches in a different folder

files1 <- list.files("index4files", full.names=T)
files1 <- files1[grep(".fastq", files1)]

files4 <- list.files("index1files", full.names=T)
files4 <- files4[grep(".fastq", files4)]

overlappingones <- read.table("index14overlapfile", stringsAsFactors=T)

for(i in 1:length(overlappingones[,1])){
	index1 <- files1[grep(overlappingones[i,1], files1)]
	index4 <- files4[grep(overlappingones[i,1], files4)]
	system(paste("cat ", index1, " >> ", index4, sep=""))
}

###Rename a bunch of files without deleting the old ones
### First make a directory for where the new files will go, then move to the directory where the old files are. Then make a file that has the column of new names next to a column of old adapters
setwd("/global/home/kebi/data/Alex/startingover_8Oct2015/IndContig")
files <- list.files("/global/home/kebi/data/Alex/startingover_8Oct2015/IndContig")
files <- files[grep(".fastq", files)]
conversion <- read.table("/global/home/kebi/data/Alex/startingover_8Oct2015/IndContig/conversion", sep = "\t", stringsAsFactors = T)


for(i in 1:nrow(conversion)){
	newfile <- files[grep(conversion[i,2],files)]
	newname <- sub(pattern=conversion[i,2], replace = conversion[i,1], newfile)
	system(paste("cp ", files[grep(conversion[i,2],files)] , " ./newfilenames/", newname, sep=""))
	
}

### Move certain files in a list to some directory
files <- list.files("directorywherethefilesaremovingfrom")
files <- files[grep("*.fastq", files)]
movethesefiles <- read.table("listoffilestomove", stringsAsFactor=T)

for(i in 1:nrow(movethesefiles)){
system(paste("cp ", files[grep(movethesefiles[i,], files)], " /global/home/kebi/data/Alex/startingover_8Oct2015/IndContig", sep=""))
}


### use ANGSD to convert a bunch of BAMs to FASTA
files <- list.files("~/data/Alex/startingover_8Oct2015/IndContig/Crotaphytus/individual_contigs/PopGen_reference/alignment")
files <- files[grep("*.bam", files)]
files <- files[seq(1,92,by=2)]
names<-as.vector(sapply(files, function(x) unlist(strsplit(x, split="_"))[2]))

for (i in 1:length(files)){
	system(paste("angsd ", "-i ", files[i], " -doFasta 1 ", "-out ", names[i], sep=""))
}

