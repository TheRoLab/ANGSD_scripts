#Wrappers for calculating Fst from multiple populations using ANGSD. To be run in R through bash wherever you are running ANGSD
#Warning: the names get a bit messy with this file, but just pay attention to the end of the filename

#start in directory with all the bamlists (alignment_out)
## Make per pop saf 
files <- list.files(getwd())
bamlist <-  files[grep("*_bamlist$", files)]

for(i in 1:length(bamlist)){
  system(paste("/global/home/kebi/programs/angsd0.913/angsd -anc ../../../../../genomes/kmer60-min500-scaffolds.fa -sites /global/home/kebi/data/Alex/voleseq/round2/demultiplexed_flipped_noclones2/alignmentout/ANGSD/keep_90percent.bed -only_proper_pairs 1 -minMapQ 2 -minQ 20 -dosaf 1 -gl 1 -bam ", bamlist[i], " -out ", bamlist[i], "_sfs", sep=""))
}

#calculate 2dsfs prior

files <- list.files(getwd())
idxlist <-  files[grep("*.saf.idx", files)]
index <- combn(8,2) #All the pairwise combinations of the 8 populations that I have
for(i in 1:ncol(index)){
  system(paste("/global/home/kebi/programs/angsd0.911/misc/realSFS ", idxlist[index[1,i]], " ", idxlist[index[2,i]]," > ",idxlist[index[1,i]],"_", idxlist[index[2,i]],".ml",sep=""))
}

## prepare the fst for window analysis

files <- list.files(getwd())
idxlist <-  files[grep("*.saf.idx$", files)]
index <- combn(8,2)
for(i in 1:ncol(index)){
  system(paste("/global/home/kebi/programs/angsd0.911/misc/realSFS fst index ", idxlist[index[1,i]], " ", idxlist[index[2,i]]," -sfs ",idxlist[index[1,i]],"_", idxlist[index[2,i]],".ml -fstout ",  idxlist[index[1,i]],"_", idxlist[index[2,i]], sep=""))
}


## get the global fst estimate. this goes into bash, not R

for pop in `ls *.fst.idx`
do
/global/home/kebi/programs/angsd0.911/misc/realSFS fst stats $pop
done

