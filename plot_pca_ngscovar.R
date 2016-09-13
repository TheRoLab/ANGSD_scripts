##plot pca from tecopa and non tecopa pops.
#this script takes the output covariance matrix from ngsCovar and makes a PCA, without blindly using the plot script from ngsCovar. 

setwd("~/Documents/Berkeley/VoleProject/bestrad2016/Analysis/aligned_radtools_noclones/")

#For a subset of my individuals
covar <- read.table("pca_tecopaonly.covar", stringsAsFactors = FALSE) 
eig <- eigen(covar, symm=TRUE)
eig$val <- eig$val/sum(eig$val)
cat(signif(eig$val, digits=3)*100,"\n") #this is the amount of variation that each axis explains

PC <- as.data.frame(eig$vectors)
PC$pop <- as.factor(c(rep("10",3),rep("11",2),rep("1",9),rep("21",4), rep("39", 7),rep("55",4), rep("9", 6), rep("F", 10), rep("Col",11)))
PC$color <- c(rep("green",3),rep("red",2),rep("orange",9),rep("blue",4), rep("pink", 7),rep("lightblue",4), rep("purple", 6), rep("black", 10), rep("grey",11))

#Plot by population
plot(PC$V1~PC$V2, col=PC$color, pch=19, ylab="PC1 (7.67%)", xlab ="PC2 (3.22%)")
legend("topleft", col=unique(PC$color), legend=c("10","11","1","21","39","55","9","Colony","F2+'F3' Colony"), pch=19, cex=0.85)




##PCA for all inds
covar <- read.table("pca_allinds.covar", stringsAsFactors = F) 
dim(covar)
eig <- eigen(covar, symm=TRUE)
eig$val <- eig$val/sum(eig$val)
cat(signif(eig$val, digits=3)*100,"\n")


PC <- as.data.frame(eig$vectors)
PC$pop <- as.factor(c(rep("10",3),rep("11",2),rep("1",9),rep("21",4), rep("39", 7),rep("55",4), rep("9", 6), rep("Col", 11), rep("DS",7), rep("F",10), rep("LE", 4), rep("LI", 4), rep("LS",5), rep("Vic",6)))
PC$color <- c(rep("green",3),rep("red",2),rep("orange",9),rep("blue",4), rep("pink", 7),rep("lightblue",4), rep("purple", 6), rep("black", 11), rep("grey",7),rep("aquamarine",10),rep("brown",4), rep("maroon",4), rep("chartreuse",5), rep("cyan",6))

plot(PC$V3~PC$V2, col=PC$color, pch=19, ylab="PC2", xlab ="PC1")
legend("topleft", col=unique(PC$color), legend=unique(PC$pop), pch=19, cex=0.85)






