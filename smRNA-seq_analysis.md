smRNA-seq analysis
===

# Required softwares

* bwa (v0.7.15)
* bedtools (v2.25.0)
* R (v3.5.0)
* fastq-dump (v2.9.2)

# Data

Data available as fastq files in NCBI BioProject [PRJNA427142](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA427142)

All libraries contain single end reads.

Download fastq files
```
fastq-dump --split-spot SRR6456282
fastq-dump --split-spot SRR6456283
fastq-dump --split-spot SRR6456284
fastq-dump --split-spot SRR6456285
fastq-dump --split-spot SRR8069222
fastq-dump --split-spot SRR8069223
fastq-dump --split-spot SRR8069224
```

Rename fastq files

```
mv SRR6456282.fastq 1794_C.fastq
mv SRR6456283.fastq 1794_D.fastq
mv SRR6456284.fastq	1794_A.fastq
mv SRR6456285.fastq	1794_B.fastq
mv SRR8069222.fastq	3634_C.fastq
mv SRR8069223.fastq	3634_B.fastq
mv SRR8069224.fastq	3634_A.fastq
```

Library codes:

* 1794_A: Block C IR #15-2
* 1794_B: Block C IR #15-3
* 1794_C: Block C IR #27-3
* 1794_D: Block C IR #27-4
* 3634_A: Col-0
* 3634_B: Block E IR #16-5
* 3634_C: Block E IR #18-5

# Read mapping

Download TAIR10 fasta reference file [here](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas).

Index the fasta file for bwa mapping
```
bwa index -p TAIR10_chr_all_bwaidx -a bwtsw TAIR10_chr_all.fas
```

Map to TAIR10 the reads of the four fastq files (4 libraries)
```
for i in *.fastq; do
    bwa aln -t 4 TAIR10_bwaidx $i | bwa samse TAIR10_bwaidx - $i | samtools view -bS - | samtools sort - $i
done
```

# Generate density plot (Fig. 2a, Extended Fig. 6b)

## Extract reads mapping at IR

Convert all bam files into bed files
```
for i in *bam; do
    bedtools bamtobed -i $i > $i.bed
done
```

Create bed files with the coordinates of the IR target +/- 100 bp to get the putative flanking reads

Block C IR

```
echo -e “chr1\t24325688\t24326434\ttarget_BlockC_plus_minus_100bp” > target_BlockC_plus_minus_100bp.bed
```

Block E IR

```
echo -e "chr1\t24334909\t24335722\ttarget_BlockE_plus_minus_100bp" target_BlockE_plus_minus_100bp.bed
```

Use bedtools to retrieve all reads overlapping the coordinates of the IR target +/- 100 bp


For Block C IR lines
```
for i in 1794*bam.bed; do
	bedtools intersect -wa -a $i -b target_BlockC_plus_minus_100bp.bed > ${i}.target_region_plus_minus_100bp.bed
done
```

For Block E IR lines and WT
```
for i in 3634*bam.bed; do
	bedtools intersect -wa -a $i -b target_BlockE_plus_minus_100bp.bed > ${i}.target_region_plus_minus_100bp.bed
done
```

Count number of reads for each libraries
```
wc –l *target_ Block*_region_plus_minus_100bp.bed
```

Only 2-6 reads are detected in the two libraries corresponding to the non-transgenic siblings (B = #15-3, C = #27-3)

To normalize read count at IR to the total number of reads per library, The number of reads in \*target_region_plus_minus_100bp.bed is divided the by number of reads in bed files for the whole library

```
for i in *fastq.bam.bed
do
	wc –l $i
done
```
```
9949902 1794_A.fastq.bam.bed
10366797 1794_B.fastq.bam.bed
10124258 1794_C.fastq.bam.bed
9998549 1794_D.fastq.bam.bed
5301603 3634_A.fastq.gz.bam.bed
5168771 3634_B.fastq.gz.bam.bed
4885825 3634_C.fastq.gz.bam.bed
```

All libraries have about 10M reads, use the values of read count for library A and D in R script for figure 2a (see below).

Bed files ready to be analyzed in R (see script below)

## Generate density plot in R (figure 2a and extended figure 6b)

### For libraries project 1794 (Block C IR)

Load R libraries

```{r}
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(Biostrings)
```

Load bed files for the 2 libraries which map smRNAs at Block C (A, D)

```{r}
readsA <- import.bed(con="1794_A.fastq.bam.bed.target_region_plus_minus_100bp.bed")
readsD <- import.bed(con="1794_D.fastq.bam.bed.target_region_plus_minus_100bp.bed")
```

Calculate coverage
```{r}
covA <- coverage(readsA)
covD <- coverage(readsD)
```

#### Plot read density for Block C (figure 2a)

Plot read density at the target region. Get first and last position for each bed file looking at 5' coordinate of first read and 3' coordinate of last read (avoid error while plotting) 
To normalize the amount of reads, divide coverage by the total number of reads of each library

```{r}
par(mfrow=c(1,1))
plot(covD[["chr1"]][24325772:24326384]/ 9998549, type="n", main="Read distribution at the target region", xlab="position",ylab="read density")
lines(covA[["chr1"]][24325772:24326336]/ 9949902, lty=2, lwd=2)
lines(covD[["chr1"]][24325772:24326384]/ 9998549,  lty=1, lwd=2)
legend(0, 6500, c("#27-4","#15-2"),lwd=c(2,2), lty=c(1,2))
```

Knowing Block C IR is located chr1:24325788..24326334, I can identify the beginning and the end of the IR on the plot

start = 24325788 - 24325772 = 16

Total length plot is 24326384 - 24325772 = 612

end = 24326384 - 24326334 = 50 so 612 - 50 = 562

Create vertical lines at the beginning and the end of the IR

```{r}
abline(v=16, col="red")
abline(v=562, col="red")
```

Position of Block C, chr1::24325922..24326268  

start = 24325922 - 24325772 = 150

end = 24326384 - 24326268  = 116 so 612 - 116 = 496

Create vertical lines at the beginning and the end of Block C

```{r}
abline(v=150, col="blue")
abline(v=496, col="blue")
```


### For libraries project 3634 (Block E IR + WT)

Load bed files for the WT Col-0 (A) and the 2 libraries which map smRNAs at Block E (B, C)

```{r}
readsA <- import.bed(con="3634_A.fastq.gz.bam.bed.target_region_BlockE_plus_minus_100bp.bed")
readsB <- import.bed(con="3634_B.fastq.gz.bam.bed.target_region_BlockE_plus_minus_100bp.bed")
readsC <- import.bed(con="3634_C.fastq.gz.bam.bed.target_region_BlockE_plus_minus_100bp.bed")
covA <- coverage(readsA)
covB <- coverage(readsB)
covC <- coverage(readsC)
```

Calculate coverage

```{r}
covA <- coverage(readsA)
covD <- coverage(readsD)
```

#### Plot read density for Block E (extended figure 6b)

Plot read density at the target region. Get first and last position for each bed file looking at 5' coordinate of first read and 3' coordinate of last read (avoid error while plotting) 
To normalize the amount of reads, divide coverage by the total number of reads of each library

```{r}
par(mfrow=c(1,1))
options("scipen"=-1) # Set y axis with scientific annotation (exponential)
plot(covC[["chr1"]][24335009:24335620]/4885825, type="n", main="Read distribution at the target region", xlab="Position (bp)",ylab="Normalized read density",  ylim=c(0,3.5e-4))
lines(covA[["chr1"]][24335009:24335674]/5301603, lty=2, lwd=2)
lines(covB[["chr1"]][24335008:24335620]/5168771, lty=1, lwd=2)
lines(covC[["chr1"]][24335009:24335620]/4885825,  lty=3, lwd=2)
legend(25, 0.00035, c("#Col-0","#16-5","#18-5"),lwd=c(2,2), lty=c(2,1,3))
```

Position of BlockE_IR = chr1:24335009-24335622 

start = 24335009 - 24335008 = 1

end =  24335622 - 24335008  = 614

```{r}
abline(v=1, col="blue")
abline(v=614, col="blue")
```

# Generate plot frequency distribution of read size (Fig. 2b, Extended Fig. 6c)

## Get read size for Block C IR

Create a report with the number of read mapping at Block C. Each line in the bam file matches with one mapped read

```
for i in *bam; do
    samtools view $i "chr1:24325922-24326301" > ${i}_reads_BlockC.sam
done
```

For the transgene region (a bit wider than the BlockC itself)

```
for i in *bam; do
    samtools view $i "chr1:24325788-24326334" > ${i}_reads_target_region.sam
done
```

Get rid of unremoved part of adapter sequences on the reads (unproper trimming of Cutadapt probably caused by short length of the reads). Considered the mapped part of each read. To do so, the following for loop takes the 19th column of a BAM file, separates the fields of the 19th column which are separated by colons, remove all '^' special characters, letter+number before and after the length of the mapped part of the reads, puts the content in *ready file. I ensure with grep that I have only values containing 2 figures. The *ready file contains the read length for each mapped read. This problem was solved for the project 3634 (Block E IR).

```
 for i in *sam;  do
    cut -f19 $i | awk -F":" '{print $3}' - | sed -e 's/\^//g' -e 's/^[0-9][A-Z][0-9][A-Z]//g' -e 's/^[0-9][A-Z]//g' -e 's/[A-Z][0-9][0-9][A-Z][0-9][0-9]$//g' -e 's/[A-Z][0-9][0-9][A-Z][0-9]$//g' -e 's/[A-Z][0-9][0-9]$//g' -e 's/[A-Z][0-9][A-Z][0-9]$//g' -e 's/[A-Z][0-9]$//g' - | grep -w '[0-9][0-9]' - > ${i}.ready
 done
```
 
*ready files can be analysed in R (see script below)*


### Plot read size distribution (figure 2b)

For normalized read length distribution, just give percent of smRNA-seq in each category as we look only into the reads into BlockC. To do, use prop.table function of contigency table generated with table()

```{r}
lengthA <- scan("1794_A.fastq.bam_reads_target_region.ready")
lengthB <- scan("1794_B.fastq.bam_reads_target_region.ready")
lengthC <- scan("1794_C.fastq.bam_reads_target_region.ready")
lengthD <- scan("1794_D.fastq.bam_reads_target_region.ready")

dataList <- list("#15-2"=lengthA, "#15-3"=lengthB, "#27-3"=lengthC, "#27-4"=lengthD)

par(mfrow=c(1,2),oma = c(0, 0, 2, 0))

plot(prop.table(table(dataList[[1]])),xlim=c(18,25),  xlab="size (nt)", ylab="number of reads", main = "transgenic #15-2", cex.main=1)
plot(prop.table(table(dataList[[4]])), xlim=c(18,25),  xlab="size (nt)", ylab="number of reads", main = "transgenic #27-4", cex.main=1)
mtext("Distribution of read sizes at target region", outer = TRUE, cex = 1)

```


## Get read size for Block E IR

Create a report with the number of read mapping at Block E. Each line in the bam file matches with one mapped read

```
for i in *bam; do
    samtools index -b $i
	samtools view $i "chr1:24335009-24335622" > ${i}_reads_BlockE.sam # to concatenate a variable, put the variable name between curly brackets
done
```

No problem of unremoved adapter sequences for the libraries of Block E IR and Col-0 but still do this step to get only the part of the reads that match the reference sequence.

```
 for i in *sam;  do
    cut -f19 $i | awk -F":" '{print $3}' - | sed -e 's/\^//g' -e 's/^[0-9][A-Z][0-9][A-Z]//g' -e 's/^[0-9][A-Z]//g' -e 's/[A-Z][0-9][0-9][A-Z][0-9][0-9]$//g' -e 's/[A-Z][0-9][0-9][A-Z][0-9]$//g' -e 's/[A-Z][0-9][0-9]$//g' -e 's/[A-Z][0-9][A-Z][0-9]$//g' -e 's/[A-Z][0-9]$//g' - | grep -w '[0-9][0-9]' - > ${i}.ready
 done
```

*ready files can be analysed in R (see script below)*

### Plot read size distribution (extended figure 6c)

For normalized read length distribution, just give percent of smRNA-seq in each category as we look only into the reads into BlockC. To do, use prop.table function of contigency table generated with table()

```{r}

lengthB <- scan("3634_B.fastq.gz.bam_reads_BlockE.sam.ready")
lengthC <- scan("3634_C.fastq.gz.bam_reads_BlockE.sam.ready")

dataList <- list("#16-5"=lengthB, "#18-5"=lengthC)

par(mfrow=c(1,2),oma = c(0, 0, 2, 0))

plot(prop.table(table(dataList[[1]])),xlim=c(18,25), ylim=c(0,0.65),  xlab="size (nt)", ylab="", main = "transgenic #16-5", cex.axis=1.3, cex.main=1)
plot(prop.table(table(dataList[[2]])), xlim=c(18,25),  ylim=c(0,0.65), xlab="size (nt)", ylab="", main = "transgenic #18-5", cex.axis=1.3, cex.main=1)
mtext("Distribution of read sizes at target region", outer = TRUE, cex = 1)

```



