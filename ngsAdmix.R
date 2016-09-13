#Instructions for doing NGSadmix with ANGSD

#Create beagle file from your bamlist:
angsd -bam tecopa_bamlist_reduced -out ngsadmix_tecopa_only -sites ../ANGSD/keep_90percent.bed -anc ../../../../../genomes/kmer60-min500-scaffolds.fa -ref ../../../../../genomes/kmer60-min500-scaffolds.fa -only_proper_pairs 1 -minMapQ 2 -minQ 20 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -doGlf 2 -nThreads 8 -GL 1 

#Loop within a loop to run each K value 9 times. (For Evanno 2005)
for(j in 1:9){
for(i  in 1:9){
  system(paste("NGSadmix -likes ngsadmix_tecopa_only.beagle.gz -P 8 -minMaf 0.05 -K ", i, " -o ", "tecopa_vicK", i,"run",j, " -minInd 38", sep=""))
}}

#catenate all of the ML values of all of the runs into a log file for submission to CLUMPAK for Evanno 2005 calculations. In terminal do:

(for log in `ls *.log`; do grep -Po 'like=\K[^ ]+' $log; done) > logfile

#Import the logfile to R, add a column to the left of each likelihood score for which K it represented
logs <- as.data.frame(read.table("logfile"))
logs$K <- c(rep("1",9),rep("2",9),rep("3",9),rep("4",9),rep("5",9),rep("6",9),rep("7",9),rep("8",9),rep("9",9))
write.table(logfile,"filelocation", row.names=F,col.names=F,quote=F)





#plot ngsAdmix outputs. You can also use CLUMPAK to do this.

admix<- t(as.matrix(read.table("~/Documents/Berkeley/VoleProject/bestrad2016/Analysis/aligned_radtools_noclones/k90_outputs/ngsAdmix/tecopa_onlyK2run1.qopt")))
dim(admix)
pops <- c(rep("10",3,), rep("11",2), rep("1",8), rep("21",4), rep("39",7), rep("55",4), rep("9",3), rep("F",10), rep("Col",2))

barplot(admix,col=1:2,space=0,border=NA,ylab="Admixture",xlab="Populations", main="Reduced Dataset Voles (K=2)")
text(c(1,4,8.5,15,20,26,30,35,42), -0.05,unique(pops),xpd=T)
abline(v=c(3,5,13,17,24,28,31,41), lty=5, lwd=2, col="white")

#Change the order manually to put the populations with the same color next to each other:
geolocated <- admix[,c(4:5,18:24,25:28,29:31,14:17,1:3,6:13,32:41,42:43)]
barplot(geolocated,col=1:2,space=0,border=NA,ylab="Admixture",xlab="Populations", main="Reduced Dataset Voles (K=2)")
abline(v=c(2,9,13,16,20,23,31,41), lty=5, lwd=2, col="white")
pops <- c(rep("11",2), rep("39",7), rep("55",4),rep("9",3), rep("21",4),rep("10",3,), rep("1",8),     rep("F",10), rep("Col",2))
text(c(1,6,11,14.5,18,22,28,36,42), -0.05,unique(pops),xpd=T)
