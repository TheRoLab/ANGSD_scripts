#Import output of ngsrelate for all individuals
setwd("~/Documents/Berkeley/VoleProject/bestrad2016/Analysis/aligned_radtools_noclones/k90_outputs/relatedness/")
relate <- read.table("relatedness_allinds_freqflipped", sep="\t",header=T)
dim(relate) #1953 pairwise comparisons -- 1953*2 = 3906

#make one column for individual 1 and one column for individual 2 in the pairwise comparisons
ind_comparisons <- unlist(strsplit(as.character(relate$Pair), split=":"))[c(seq(from=1, to=3906,by=2))] #extract the individuals that were compared
ind <- unlist(strsplit(as.character(ind_comparisons), split=","))
no_parentheses <- gsub('\\(',"", ind)
no_parentheses <- gsub('\\)',"", no_parentheses)
ind1 <- no_parentheses[c(seq(from=1, to=3905,by=2))]
ind2 <- no_parentheses[c(seq(from=2, to=3906,by=2))]
#make the inds 1 based instead of 0 based. Now individual 1 is the first individual in your bamlist
ind1 <- as.numeric(ind1)+1
ind2 <- as.numeric(ind2)+1
relate_withinds <- data.frame(ind1,ind2,relate$k0,relate$k1,relate$k2) 



#Mean of all K0, K1 and K2 where individual from Marsh 10 was compared to an individual from Marsh 10. Individuals from Marsh 10 are individuals 1, 2 and 3 in the bamlist

marshten <- relate_withinds[which(relate_withinds$ind1 == 1 | relate_withinds$ind1 == 2 | relate_withinds$ind1 == 3),]
marshten <- marshten[which(marshten$ind2 == 1 | marshten$ind2 == 2 | marshten$ind2 == 3),]

#Co-ancestry coefficient, theta = k1/4 + k2/2. The mean is:
 mean((marshten$relate.k1/4) + (marshten$relate.k2/2))
 0.1048963
 #repeated for each pop:
 
eleven <- relate_withinds[which(relate_withinds$ind1 == 4 | relate_withinds$ind1 == 5),]
eleven <- eleven[which(eleven$ind2 == 4 | eleven$ind2 == 5),]
mean((eleven$relate.k1/4) + (eleven$relate.k2/2))
0.291638

one <- relate_withinds[which(relate_withinds$ind1 == 6 | relate_withinds$ind1 == 7 | relate_withinds$ind1 == 8 | relate_withinds$ind1 == 9 | relate_withinds$ind1 == 10 | relate_withinds$ind1 == 11 |relate_withinds$ind1 == 12 | relate_withinds$ind1 == 13),]
one <- one[which(one$ind2 == 6 | one$ind2 == 7 | one$ind2 == 8 | one$ind2 == 9 | one$ind2 == 10 | one$ind2 == 11 |one$ind2 == 12 | one$ind2 == 13),]
one
mean((one$relate.k1/4) + (one$relate.k2/2))
0.08324853

twentyone <- relate_withinds[which(relate_withinds$ind1 == 14 | relate_withinds$ind1 == 15 | relate_withinds$ind1 == 16 | relate_withinds$ind1 == 17),]
twentyone <- twentyone[which(twentyone$ind2 == 14 | twentyone$ind2 == 15 | twentyone$ind2 == 16 | twentyone$ind2 == 17),]
mean((twentyone$relate.k1/4) + (twentyone$relate.k2/2))
0.1924835



thirtynine <- relate_withinds[which(relate_withinds$ind1 == 18 | relate_withinds$ind1 == 19 | relate_withinds$ind1 == 20 | relate_withinds$ind1 == 21 | relate_withinds$ind1 == 22 | relate_withinds$ind1 == 23 |relate_withinds$ind1 == 24),]
thirtynine <- thirtynine[which(thirtynine$ind2 == 18 | thirtynine$ind2 == 19 | thirtynine$ind2 == 20 | thirtynine$ind2 == 21 | thirtynine$ind2 == 22 | thirtynine$ind2 == 23 |thirtynine$ind2 == 24),]
mean((thirtynine$relate.k1/4) + (thirtynine$relate.k2/2))
0.1138538


fiftyfive <- relate_withinds[which(relate_withinds$ind1 == 25 | relate_withinds$ind1 == 26 | relate_withinds$ind1 == 27 | relate_withinds$ind1 == 28),] 
fiftyfive <- fiftyfive[which(fiftyfive$ind2 == 25 | fiftyfive$ind2 == 26 | fiftyfive$ind2 == 27 | fiftyfive$ind2 == 28),] 
mean((fiftyfive$relate.k1/4) + (fiftyfive$relate.k2/2))
0.2568861

nine <- relate_withinds[which(relate_withinds$ind1 == 29 | relate_withinds$ind1 == 30 | relate_withinds$ind1 == 31),] 
nine <- nine[which(nine$ind2 == 29 | nine$ind2 == 30 | nine$ind2 == 31),] 
mean((nine$relate.k1/4) + (nine$relate.k2/2))
0.1454075

F <- relate_withinds[which(relate_withinds$ind1 == 32 |relate_withinds$ind1 == 33 |relate_withinds$ind1 == 34 | relate_withinds$ind1 == 35 | relate_withinds$ind1 == 36 | relate_withinds$ind1 == 37 | relate_withinds$ind1 == 38 | relate_withinds$ind1 == 39 |relate_withinds$ind1 == 40 | relate_withinds$ind1 == 41),]
F <- F[which(F$ind2 == 32 |F$ind2 == 33 |F$ind2 == 34 | F$ind2 == 35 | F$ind2 == 36 | F$ind2 == 37 | F$ind2 == 38 | F$ind2 == 39 |F$ind2 == 40 | F$ind2 == 41),]
mean((F$relate.k1/4) + (F$relate.k2/2))
0.12013


Col <- relate_withinds[which(relate_withinds$ind1 == 42 |relate_withinds$ind1 == 43 |relate_withinds$ind1 == 44 |relate_withinds$ind1 == 45 | relate_withinds$ind1 == 46 | relate_withinds$ind1 == 47 |relate_withinds$ind1 == 48 |relate_withinds$ind1 == 49 |relate_withinds$ind1 == 50 |relate_withinds$ind1 == 51 |relate_withinds$ind1 == 52 |relate_withinds$ind1 == 53),]
Col <- Col[which(Col$ind2 == 42 |Col$ind2 == 43 |Col$ind2 == 44 |Col$ind2 == 45 | Col$ind2 == 46 | Col$ind2 == 47 |Col$ind2 == 48 |Col$ind2 == 49 |Col$ind2 == 50 |Col$ind2 == 51 |Col$ind2 == 52 |Col$ind2 == 53),]
mean((Col$relate.k1/4) + (Col$relate.k2/2))
0.09772755


#For all of tecopa (plus a few from the Colony)
mean((relate_withinds$relate.k1/4)+(relate_withinds$relate.k2/2))
0.058369 #doesn't make much sense
#mean of all of the above values, excluding the colony, is 0.1866

####Try again for Pop 9 and Vic (should be decidedly less related because they are hundreds of miles apart)

relate <- read.table("nine_vic_relate.res",header=T)
dim(relate) #66 pairwise comparisons
ind_comparisons <- unlist(strsplit(as.character(relate$Pair), split=":"))[c(seq(from=1, to=132,by=2))] #extract the individuals that were compared
ind <- unlist(strsplit(as.character(ind_comparisons), split=","))
no_parentheses <- gsub('\\(',"", ind)
no_parentheses <- gsub('\\)',"", no_parentheses)
ind1 <- no_parentheses[c(seq(from=1, to=131,by=2))]
ind2 <- no_parentheses[c(seq(from=2, to=132,by=2))]
#make the inds 1 based instead of 0 based. Now individual 1 is the first individual in your bamlist
ind1 <- as.numeric(ind1)+1
ind2 <- as.numeric(ind2)+1
relate_withinds <- data.frame(ind1,ind2,relate$k0,relate$k1,relate$k2) 

#relatedness should be 0 because they're hundreds of miles apart
nine <- relate_withinds[which(relate_withinds$ind1 == 1 | relate_withinds$ind1 == 2 | relate_withinds$ind1 == 3 | relate_withinds$ind1 == 4 | relate_withinds$ind1 == 5 | relate_withinds$ind1 == 6),]
nine <- nine[which(nine$ind2 == 7 | nine$ind2 == 8 | nine$ind2 == 9 | nine$ind2 == 10 | nine$ind2 == 11),]
mean((nine$relate.k1/4) + (nine$relate.k2/2))
0 


#What about within Victorville? 

vic <- relate_withinds[which(relate_withinds$ind1 == 7 | relate_withinds$ind1 == 8 | relate_withinds$ind1 == 9 | relate_withinds$ind1 == 10 | relate_withinds$ind1 == 11),]
vic <- vic[which(vic$ind2 == 7 | vic$ind2 == 8 | vic$ind2 == 9 | vic$ind2 == 10 | vic$ind2 == 11),]

mean((vic$relate.k1/4) + (vic$relate.k2/2))
0.09399333

#looks good!
