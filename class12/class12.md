---
title: "class12"
author: "Michelle Janossy"
date: "5/9/2019"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Setup for Docking 

We will first prepare our HIV-Pr system for drug docking by making a protein only PDB format file (i.e. we will remove water, existing ligand, etc)

```{r}
library(bio3d)
file <- get.pdb("1hsg")
```

```{r}
pdb <- read.pdb(file)
pdb
```

```{r}
prot <- atom.select(pdb, "protein", value = TRUE)
write.pdb(prot, file = "1hsg_protein.pdb")
prot
```

```{r}
lig <- atom.select(pdb, "ligand", value = TRUE)
write.pdb(lig, file = "1hsg_ligand.pdb")
lig
```

We will load these into ADT to add hydrogens and make PDBQT files for each 

Now we can run autodoc vina with the cmd

> "\Program Files (x86)\The Scripps Research Institute\Vina\vina.exe" --
config config.txt --log log.txt

```{r}
library(bio3d)
res <- read.pdb("all.pdbqt", multi = TRUE)
write.pdb(res, "results.pdb")
```

```{r}
#res <- read.pdb("all.pdbqt", multi = TRUE)
ori <- read.pdb("1hsg_ligand.pdbqt")
rmsd(ori,res)
```

## Normal mode analysis for flexibility prediction

```{r}
pdb <- read.pdb("1hel")
modes <- nma( pdb )
m7 <- mktrj(modes, mode=7, file="mode_7.pdb")
```













