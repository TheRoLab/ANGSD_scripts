
import pandas as pd
import numpy as np

file=input("Please give the full path to the .pestPG output file: ")

df = pd.DataFrame.from_csv(file, sep="\t")

watterson = pd.DataFrame.mean(df["tW"])
pairwise = pd.DataFrame.mean(df["tP"])
Tajima = pd.DataFrame.mean(df["Tajima"])

Dattempt = (pairwise - watterson)/np.sqrt(pd.DataFrame.var(df["tP"]-df["tW"]))

print("Watterson's Theta =  ", watterson)
print("Pairwise Pi a.k.a. Theta pi =  ", pairwise)
print("Mean Tajima's D using Alex's Method =  ", Dattempt)
print("Mean Tajima's D using ANGSD Method =  ", Tajima)

