## shRNA counting for pooled shrna-seq screening
### R code from vignette source 'pooledScreenAnalysis.rnw'
```
setwd("~/shRNA-seq/MSC_shRNA_021116/exp1/edgeR_method/")
getwd()
rm(list=ls())
options(digits=3)
options(width=90)
library(edgeR)
```

Read in sample & hairpin information
```
sampleanno = read.table("samples1.txt", sep=",", header=TRUE)
sampleanno

hairpinseqs = read.table("library_blankremoved_depulicated.txt", sep="\t", header=TRUE)
hairpinseqs[,3]
```

To convert columns of ID and sequences from factors to character
```
class(hairpinseqs$ID)
hairpinseqs[] <- lapply(hairpinseqs, as.character)
nchar(hairpinseqs$Sequences)
```

have to remove all (shRNA) seqences that are not 19 bp
```
data = hairpinseqs[(which(nchar(hairpinseqs$Sequences) == 19)),]
data[1:5,]
write.table(data, "data.txt", sep = "\t")
```

Process raw sequences from fastq file
```
x = processAmplicons("all_reads_R2.fastq", barcodefile="samples1.txt", 
                     hairpinfile="data.txt", hairpinStart=58, hairpinEnd=76,allowShifting=TRUE, shiftingBase=5,verbose=TRUE)
```


Filter hairpins with low counts, countspermillion > 0.5, and rowSums > 3
```
sel = rowSums(cpm(x$counts)>0.5)>=3
x = x[sel,]
x
```
Plot number of hairpins that could be matched per sample
```
par(mfrow=c(2,1))
barplot(colSums(x$counts), las=2, main="Counts per index", cex.names=0.5, cex.axis=0.8, ylim=c(0,3000000))
colSums(x$counts)
```

Plot per hairpin totals across all samples
```
barplot(rowSums(x$counts), las=2, main="Counts per hairpin", cex.names=0.5, cex.axis=0.8, ylim=c(0,400000))
barplot(rowSums(x$counts), las=2, main="Counts per hairpin", cex.names=0.5, cex.axis=0.8, ylim=c(0,200000))
```
subset the counts, genes and other info from the x dataframe
```
xcounts <- (x$counts)
xgenes <- (x$genes)
x$samples$group
lib_size = as.vector(x$samples$lib.size)
name_groups = as.vector(x$samples$group)
names(lib_size) <- name_groups
```

look at the first 5 lines
```
xcounts[1:5,]
xcounts2=xcounts
```

rename column
```
colnames(xcounts) <- c("baseline", "p16", "p17")
```

extract all hairpin IDs from count table
```
hairpinIDs <- rownames(xcounts2)
length(hairpinIDs)
```
split each row of the xcount2 table into lists
```
sample_counts <- split(xcounts2, row(xcounts2))
```
add names of hairpinIDs to list
```
names(sample_counts)<- hairpinIDs
sample_counts[1:5]
```
output the gene count, gene ID and library size file
```
write.csv(xcounts, file = "counts.csv")
write.csv(xgenes, file = "genes.csv")
write.csv(lib_size, file = "lib_size.csv")
```
