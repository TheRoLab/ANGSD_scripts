##Check for missing data from -doGeno 2 output.

gen.mat<- read.table("~/Documents/Berkeley/VoleProject/bestrad2016/Analysis/aligned_radtools_noclones/calledgenos_min_maxdepth.geno")
head(gen.mat)
gen.mat <- gen.mat[,3:84] #remove site and position columns
gen.matrix.raw <- t(gen.mat)
head(gen.matrix.raw)
dim(gen.mat) #135564 SNPs
dim(gen.matrix.raw)

samples <- gen.matrix.raw
samples[samples >= 0] <- 2 #If we measured a SNP, give it a 2
samples[samples == -1] <- 0 #If the SNP is missing, give it a 0
head(samples)
dim(samples)
hist(rowSums(samples)/135564) #Histogram of percent missing data
sort(rowSums(samples)/135564) #shows the individuals with the most missing data (corresponding to the order of individuals in your bamlist)


#Make a PCA from called SNPs (a la Gideon Bradburg)

#First, subset out data with more than 50% missing data:
gen.mat.sub <- gen.matrix.raw[c(1:8,10:29,31,34:46,48:63,65:66,68:82),]
subsetted <-gen.mat.sub

#Function to make SNPs random with respect to whether the homozygous is ancestral or derived
random.switcharoo <- function(x,sample.size=2){
  #recover()
  x <- ifelse(rep(runif(1) < 0.5,length(x)),
              x,
              sample.size-x)
  return(x)
}
switcharoo.data <- function(frequencies){
  frequencies <- apply(frequencies,2,random.switcharoo)
  return(frequencies)
}
gen.matrix.switcharoo <- switcharoo.data(gen.mat.sub)
gen.matrix.switcharoo[which(subsetted==-1)] <- 0

samples <- subsetted
samples[samples >= 0] <- 2
samples[samples == -1] <- 0

hist(rowSums(samples)/135564)

#PCA based on allele frequencies
sample.freqs <- gen.matrix.switcharoo/samples
sample.cov <- cov(t(sample.freqs),use="pairwise.complete.obs")
eigen.obj <- eigen(sample.cov)

PCA <- as.data.frame(eigen.obj$vectors)

PCA$pop <- as.factor(c(rep("10",3),rep("11",2),rep("1",8),rep("21",4), rep("39", 7),rep("55",4), rep("9", 3), rep("Col", 11), rep("DS",6), rep("F",10), rep("LE", 2), rep("LI", 4), rep("LS",5), rep("Vic",6)))
PCA$color <- c(rep("green",3),rep("red",2),rep("orange",8),rep("blue",4), rep("pink", 7),rep("lightblue",4), rep("purple", 3), rep("black", 11), rep("grey",6),rep("aquamarine",10),rep("brown",2), rep("maroon",4), rep("chartreuse",5), rep("cyan",6))


plot(PCA$V1, PCA$V2, col=PCA$color, pch=19, xlab="PC1", ylab="PC2")