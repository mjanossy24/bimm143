Class 18
================
Michelle Janossy
5/30/2019

\#\#Exploring the Cancer Genome Atlas install packages we
need

``` r
BiocManager::install(c("GenomicDataCommons", "TCGAbiolinks", "maftools"))
```

load packages

``` r
#library(GenomicDataCommons)
#library(TCGAbiolinks)
#library(maftools)
```

Can we talk with the NCI-GDC?

``` r
#GenomicDataCommons::status()
```

``` r
#cases_by_project <- cases() %>%
  #facet("project.project_id") %>%
  #aggregations()
#head(cases_by_project)
```

``` r
#x <- cases_by_project$project.project_id

# Make a custom color vector for our plot
#colvec <- rep("lightblue", nrow(x))
#colvec[29] <- "red"

# Plot with 'log' for y axis and rotate labels with 'las'
#par(___)  
#barplot(x$doc_count, names.arg=x$key, log="y", col=colvec, las=2)
```

\#\#Part 2

``` r
#library(bio3d)
#seqs <- read.fasta("lecture18_sequences.fa")
#seqs
```

``` r
## Calculate positional identity scores
#ide <- conserv(seqs$ali, method="identity")
#mutant.sites <- which(ide < 1) 

## Exclude gap possitions from analysis
#gaps <- gap.inspect(seqs)
#mutant.sites <- mutant.sites[mutant.sites %in% gaps$f.inds]

#mutant.sites
```

``` r
## Make a "names" label for our output sequences (one per mutant)
#mutant.names <- paste0(seqs$ali["P53_wt",mutant.sites],mutant.sites,seqs$ali["P53_mutant",mutant.sites])

#mutant.names
```

``` r
## Sequence positions surounding each mutant site
#start.position <- mutant.sites - 8
#end.position <-  mutant.sites + 8

# Blank matrix to store sub-sequences
#store.seqs <- matrix("-", nrow=length(mutant.sites), ncol=17)
#rownames(store.seqs) <- mutant.names

## Extract each sub-sequence
#for(i in 1:length(mutant.sites)) {
  #store.seqs[i,] <- seqs$ali["P53_mutant",start.position[i]:end.position[i]]}

#store.seqs
```

``` r
## Output a FASTA file for further analysis
#write.fasta(seqs=store.seqs, ids=mutant.names, file="subsequences.fa")
```
