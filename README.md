# HyPhy-molecular-evolution---RELAX

# Preparing FASTA alignments and trees for analysis with [HyPhy](http://hyphy.org/)

## In unix, remove all stop codons from fasta files. Will add "_macse_NT" to all

````bash
macse_v1.2.jar -prog exportAlignment -align clpP-new-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rpl2-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rps12-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rpl14-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rpl16-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rpl20-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rpl23-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rpl36-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rps2-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rps3-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rps4-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rps8-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rps11-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rps14-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rps16-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rps18-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align rps19-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align ycf1-all-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN ;
macse_v1.2.jar -prog exportAlignment -align ycf2-alignment.fasta -codonForFinalStop --- -codonForInternalStop NNN
````

## Now in R, read in all fasta files and your tree to check them against each other

```{r}
library(ape)
  
clpP   <- read.dna(file = "clpP-new-alignment_macse_NT.fasta", format = "fasta", as.character = T)
rpl2   <- read.dna(file = "rpl2-alignment_macse_NT.fasta",     format = "fasta", as.character = T)
rps12  <- read.dna(file = "rps12-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
rpl14  <- read.dna(file = "rpl14-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
rpl16  <- read.dna(file = "rpl16-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
rpl20  <- read.dna(file = "rpl20-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
rpl23  <- read.dna(file = "rpl23-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
rpl36  <- read.dna(file = "rpl36-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
rps2   <- read.dna(file = "rps2-alignment_macse_NT.fasta",     format = "fasta", as.character = T)
rps3   <- read.dna(file = "rps3-alignment_macse_NT.fasta",     format = "fasta", as.character = T)
rps4   <- read.dna(file = "rps4-alignment_macse_NT.fasta",     format = "fasta", as.character = T)
rps8   <- read.dna(file = "rps8-alignment_macse_NT.fasta",     format = "fasta", as.character = T)
rps11  <- read.dna(file = "rps11-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
rps14  <- read.dna(file = "rps14-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
rps16  <- read.dna(file = "rps16-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
rps18  <- read.dna(file = "rps18-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
rps19  <- read.dna(file = "rps19-alignment_macse_NT.fasta",    format = "fasta", as.character = T)
ycf1   <- read.dna(file = "ycf1-all-alignment_macse_NT.fasta", format = "fasta", as.character = T)
ycf2   <- read.dna(file = "ycf2-alignment_macse_NT.fasta",     format = "fasta", as.character = T)

tre <- read.tree("H4-boot-hyphy.tre")

```

##  Check in R whether tree and alignment contain all identical taxa

```{r}

library(geiger)

tmp <- name.check(tre,clpP)

## If taxa are missing from an alignment, drop them from the tree
clPptree<-drop.tip(tre,c("Pecteilis_hawkesiana", "Pecteilis_susannae"))
tmp <- name.check(clPptree,clpP)
tmp

## Save the tree to use for that gene specifically
write.tree(clPtree,file="clpP.tre)
```

## 

## Now, edit trees to tag test branches
You will need to re-insert '{myco}' or '{photo}' (or whatever) after ':' following species names for all test branches in each tree!!!

{myco} = test branches, {photo} = reference branches

## If Hyphy is not installed, install it with conda

```bash
conda create -n hyphy

conda activate hyphy

conda install hyphy
```


## Run HyPhy on all datasets

```bash
conda activate hyphy   # if not activated. You will see (hyphy) and not (base) if activated

hyphy relax --alignment accD-codon_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rpl16-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rpl20-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rpl23-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rpl36-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rps2-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rps8-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rpl2-alignment_macse_NT.fasta --tree rpl2.tre
hyphy relax --alignment rps4-alignment_macse_NT.fasta --tree rps4.tre
hyphy relax --alignment rps12-alignment_macse_NT.fasta --tree rps12.tre
hyphy relax --alignment rps16-alignment_macse_NT.fasta --tree rps16.tre
hyphy relax --alignment rps18-alignment_macse_NT.fasta --tree rps18.tre
hyphy relax --alignment rpl14-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rps3-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rps19-alignment_macse_NT.fasta --tree rps19.tre
hyphy relax --alignment rps11-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rps14-alignment_macse_NT.fasta --tree H4-boot-hyphy.tre
hyphy relax --alignment rps16-alignment_macse_NT.fasta --tree rps16.tre
hyphy relax --alignment ycf1-all-alignment_macse_NT.fasta --tree ycf1.tre
hyphy relax --alignment ycf2-alignment_macse_NT.fasta --tree ycf2.tre
hyphy relax --alignment clpP-new-alignment_macse_NT.fasta --tree clpP.tre
```

## View output .json files with [HyPhy Vision](http://vision.hyphy.org/RELAX)
