## Getting help
help(rep) # Help for particular function (package must be loaded)
?rep # Help for particular function (package must be loaded)
??ANOVA # Search for the term within all installed packages
help.search("analysis of variance") # Search for the phrase within all installed packages - return list of hits sorted according to type and package (i.e. package::function)
require(sos) # More comprehensive search from packages
findFn("DAPC") # Search for function name

## Set working directory
setwd("~/dokumenty/vyuka/r_mol_data/examples/") # Create YOUR OWN empty directory and modify the path accordingly!
getwd() # Verifies where we are
dir() # Lists files and folders on the disk
ls() # Lists currently available R objects

## Get information about object classes
getClassDef("data.frame") # Or select any other class name

## Basic operations
x <- c(5, 6, 7, 8, 9) # Creates vector (see also ?rep)
x # Print "x" content
c() # Is generic function to concatenate objects into new one
length(x) # Length of the object - for matrices and DF use dim()
str(x) # Information about structure of the object
mode(x) # Gets type of storage mode of the object
class(x) # Shows class of the object
x[2] # Shows second element of the object
x <- x[-5] # Removes fifth element
x # See modified object
y <- matrix(data=5:20, nrow=4, ncol=4) # Creates a matrix
is.matrix(y) # Is it matrix? Try is.<TAB><TAB> TAB key shows available functions and objects starting by typed text
y # Prints the matrix
dim(y) # Dimensions of "y"
y[,2] # Prints second column
y[3,] # Prints third row
y[4,3] # Prints element from fourth row and third column
c(x, y[4,]) # Concatenate "x" and 4th row of "y"
x2 <- c(x, 3, 2, 1) # Concatenate "x" and values "3", "2" and "1"
x <- y[2,] # Replaces "x" by second row of "y" (no warning) - R doesn't ask neither notifies when overwriting objects! Be careful!
x # See modified object
rm(x) # Deletes x (no warning)
y[,1:3] # Prints first through third column of the matrix
y[3,] <- rep(x=20, each=4) # Replaces third line by value of 20
y # See modified object
y[y==20] <- 10 # If value of y's element is 20, replace it by 10
y # See modified object
summary(y) # Basic statistics - according to columns
colnames(y) <- c("A", "B", "C", "D") # Set column names - Objects and functions are without quotation marks; files, text with
colnames(y) # Prints column names, use rownames() in same way
y # See modified object
y[,"C"] # Prints column C (R is case sensitive!)
y[,c("C", "B")] # Extract columns "C" and "B"
t(y) # Transposes the matrix
diag(y) # Get diagonal of the matrix
diag(y) <- rep(x=c(50, 100), times=2) # Replace diagonal by repetition of values 50 and 100
y # See modified object
y <- as.data.frame(y) # Turns into DF (see other functions as.*)
class(y) # Is it data frame now?
y[y==17] <- "NA" # Removes values of 17
y # See modified object
y$B # Gets variable B of data frame y ($ works similarly in S3 objects)
# When loading saved project, you have to load again libraries and scripts (see further), data objects are restored
# This can be conveniently done in RStudio/RKWard
save(list=ls(), file="test.RData") # Saves all objects during the work
load("test.RData") # Loads saved R environment with all objects
# When loading saved project, you have to load again libraries and scripts (see further), data objects are restored
fix(y) # Use to edit matrices, data frames, functions, ...
rm(y) # Removing...

## Some basic statistics
# Basic summary statistics
summary(iris)
# Boxplot comparing sepal lengths of the three species
boxplot(formula=iris$Sepal.Length~iris$Species)
# Extract only setosa
setosa <- iris[which(iris$Species=="setosa"),]
# Testing correlation between sepal length and width of setosa
cor.test(x=setosa$Sepal.Length, y=setosa$Sepal.Width)
# Plot correlation between sepal length and width of setosa
plot(x=setosa$Sepal.Length, y=setosa$Sepal.Width, main="Sepals", xlab="Length", ylab="Width", pch=16, cex=1.5, col="red")
# Add linear model line
abline(reg=lm(formula=setosa$Sepal.Width~setosa$Sepal.Length), col='brown', lwd=3)
# Add sepal width mean and standard deviation
abline(h=mean(setosa$Sepal.Width), col="green", lwd=2, lty=2)
abline(h=mean(setosa$Sepal.Width) + sd(setosa$Sepal.Width), col="green", lwd=2, lty=3)
abline(h=mean(setosa$Sepal.Width) - sd(setosa$Sepal.Width), col="green", lwd=2, lty=3)
# ANOVA
# Basic analysis of variance
aov(formula=iris[["Sepal.Length"]]~iris[["Species"]])
# Wrapping aov() by summary() gives output which is easier to read (this is common practice for plenty of outputs)
summary(aov(formula=iris[["Sepal.Length"]]~iris[["Species"]]))

## Packages and repositories
# Set repositories
# We will need extra repositories
getOption("repos") # Shows actual repositories
options(repos=c("https://mirrors.nic.cz/R/", "https://r-forge.r-project.org/", "https://rforge.net/"))
options() # Generic function to modify various settings
?options # Gives details
# Install packages
# Installation of multiple packages may sometimes fail - install then packages in smaller groups or one by one
install.packages(pkgs=c("ade4", "adegenet", "adegraphics", "adephylo", "adespatial", "akima", "ape", "BiocManager", "caper", "corrplot", "devtools", "gee", "geiger", "ggplot2", "gplots", "hierfstat", "ips", "kdetrees", "lattice", "mapdata", "mapplots", "mapproj", "maps", "maptools", "nlme", "PBSmapping", "pegas", "permute", "phangorn", "philentropy", "phylobase", "phytools", "picante", "plotrix", "poppr", "raster", "rentrez", "rgdal", "RgoogleMaps", "Rmpi", "rworldmap", "rworldxtra", "seqinr", "shapefiles", "snow", "sos", "sp", "spdep", "splancs", "StAMPP", "TeachingDemos", "tripack", "vcfR", "vegan"), repos="https://mirrors.nic.cz/R/", dependencies="Imports")
?install.packages # See for more options
# Updates installed packages (by default from CRAN)
update.packages(ask=FALSE)
# Upgrade all packages e.g. from R 3.5 to 3.6
install.packages(pkgs=installed.packages()) # Packages installed manually from different resources must be reinstalled manually

# Install package Geneland (since version 4 not available in CRAN anymore)
# Other packages used when using Geneland
# Needed is PBSmapping or mapproj for conversion of coordinates
# GUI uses for parallelization snow and Rmpi
# RgoogleMaps (requires rgdal) can be used to plot Geneland output on top of Google map, maptools, shapefiles (requires foreign) and tripack on GIS layer
install.packages(pkgs=c("PBSmapping", "mapproj", "rgdal", "RgoogleMaps", "Rmpi", "sp", "maptools", "shapefiles", "snow", "tripack"), repos="https://mirrors.nic.cz/R/", dependencies="Imports")
# On Windows install Rtools from https://cran.r-project.org/bin/windows/Rtools/
# Install Geneland from GitHub
# Devtools package is required to install from GitHub
if(! "devtools" %in% installed.packages()) {install.packages("devtools")}
# Install Geneland itself
devtools::install_github("gilles-guillot/Geneland")

# Install package phyloch not available in any repository
# If not done already, install required packages first
install.packages(pkgs=c("ape", "colorspace", "XML"), dependencies="Imports")
# It is possible to specify direct path (local or web URL) to package source
install.packages(pkgs="http://www.christophheibl.de/phyloch_1.5-3.tar.gz", repos=NULL, type="source")

# We will load packages by library(package) one by one when needed...

## Bioconductor
# Install CRAN package BiocManager used to manage Bioconductor packages
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
?BiocManager::install # See options
# Install Bioconductor packages
BiocManager::install("muscle") # Simplest usage
BiocManager::install(pkgs=c("Biostrings", "muscle"))
# Update installed packages
BiocManager::install()
# Search for Bioconductor packages
BiocManager::available() # List everything
?BiocManager::available # See options

# Standard installation
install.packages(c("adegenet", "poppr", "phytools"))
update.packages() # Update packages

# Installation from custom repository (ParallelStructure is not used in our course)
install.packages("ParallelStructure", repos="https://r-forge.r-project.org/")
?install.packages # See help for details

## Data

# R training data
# List data available in currently loaded packages
data()
# List data available in all installed packages (can be very long)
data(package=rownames(installed.packages()))
# Load selected data set
# For example phylogeny and species traits of shorebirds from package caper (we will use it much later)
data(shorebird, package="caper")
# Optional method (load respective library and then data)
library(caper) # Library containing (among others) desired data
data(shorebird) # Loading the data
# See content of the data set
?shorebird # Or
?caper::shorebird
# R training data we will use during the course
?adegenet::microbov
?adegenet::rupica
?ape::carnivora
?caper::shorebird

## Libraries for loading data and population genetic analysis

# Load needed libraries
library(ape)
library(ade4)
library(adegenet)
library(pegas)
library(poppr)

## Work with data

# Load training data
hauss.loci <- read.loci("https://soubory.trapa.cz/rcourse/haussknechtii_ssrs.txt", header=TRUE, loci.sep="\t", allele.sep="/", col.pop=2, col.loci=3:14, row.names=1)
# Data control
hauss.loci
print(hauss.loci, details=TRUE)

# Conversion of loci to genind - used for many analysis
hauss.genind <- loci2genind(hauss.loci)
# See population names
pop(hauss.genind)
hauss.genind$pop

# Read coordinates
hauss.coord <- read.csv("https://soubory.trapa.cz/rcourse/haussknechtii_coordinates.csv", header=TRUE, sep="\t", quote="", dec=".", row.names=1)
hauss.coord

# Add coordinates - note identification of slots within object
hauss.genind$other$xy <- hauss.coord
# See result
hauss.genind$other$xy
hauss.genind

# Conversion to genpop - for population-level analysis
hauss.genpop <- genind2genpop(hauss.genind, process.other=TRUE)
# See result
hauss.genpop

# Removes missing data - see ?missingno for types of dealing with them
hauss.genind.cor <- missingno(pop=hauss.genind, type="mean", cutoff=0.1, quiet=FALSE) # Use with caution! It modifies original data!
?missingno # See other options of handling missing data
# Convert corrected genind to loci
hauss.loci.cor <- genind2loci(hauss.genind.cor)

# Writes loci file to the disk
write.loci(hauss.loci.cor, file="hauss.loci.cor.txt", loci.sep="\t", allele.sep="/")

# Read data sets from various software
?read.genalex # poppr - reads *.csv file
?read.fstat # adegenet - reads *.dat files, only haploid/diploid data
?read.genetix # adegenet - reads *.gtx files, only haploid/diploid data
?read.genepop # adegenet - reads *.gen files, only haploid/diploid data
?read.structure # adegenet - reads *.str files, only haploid/diploid data
?import2genind # adegenet - more automated version of above functions

# Import of triploid (polyploid) microsattelites
tarax3n.table <- read.table("https://soubory.trapa.cz/rcourse/tarax3n.txt", header=TRUE, sep="\t", quote="", row.names=1)
# Check the data
tarax3n.table
class(tarax3n.table)
dim(tarax3n.table)

# Convert to genind
tarax3n.genind <- df2genind(X=tarax3n.table[,1:6], sep="/", ncode=3, pop=tarax3n.table[["pop"]], ploidy=3, type="codom")
# See resulting genind object
tarax3n.genind
summary(tarax3n.genind)

# Import of AFLP data
amara.aflp <- read.table(file="https://soubory.trapa.cz/rcourse/amara_aflp.txt", header=TRUE, sep="\t", quote="")
amara.aflp
dim(amara.aflp)
class(amara.aflp) # Must be matrix or data frame
# Populations - just one column with population names for all inds
amara.pop <- read.table(file="https://soubory.trapa.cz/rcourse/amara_pop.txt", header=TRUE, sep="\t", quote="")
amara.pop
dim(amara.pop)
# Create genind object - ind.names and loc.names are taken from X
aflp.genind <- df2genind(X=amara.aflp, sep="", ind.names=NULL, loc.names=NULL, pop=amara.pop[,1], type="PA")
indNames(aflp.genind) <- amara.aflp[,1] # Add individual names
aflp.genind
summary(amara.aflp)

# Some useful functions for data manipulations
?genind2df # adegenet - export into data frame
?genind2genalex # poppr - export for genalex
?splitcombine # poppr - edits population hierarchy
?popsub # poppr - extracts only selected population(s)
?clonecorrect # poppr - corrects for clones
?informloci # poppr - removes uninformative loci
?seppop # adegenet - separates populations from genind or genlight object
?seploc # adegenet - splits genind, genpop or genlight by markers
?alleles2loci # pegas - transforms a matrix of alleles into "loci"
# seppop and seploc return lists of objects - for further analysis
# read manuals (?...) of the functions before usage

# Reading FASTA (reads also another formats, see ?read.dna), sequences of flu viruses from various years from USA (Adegenet toy data)
usflu.dna <- read.dna(file="https://adegenet.r-forge.r-project.org/files/usflu.fasta", format="fasta")
# Check the object
class(usflu.dna)
usflu.dna
# Another possibility (only for FASTA alignments, same result):
usflu.dna2 <- fasta2DNAbin(file="https://adegenet.r-forge.r-project.org/files/usflu.fasta") # Normally keeps only SNP - see ?fasta2DNAbin
# Check the object
class(usflu.dna2)
usflu.dna2
as.character(usflu.dna2)[1:5,1:10]
dim(usflu.dna2) # Does it have correct size?
# Read anotations
usflu.annot <- read.csv("https://adegenet.r-forge.r-project.org/files/usflu.annot.csv", header=TRUE, row.names=1)
head(usflu.annot)
# Convert DNAbin to genind - only polymorphic loci are retained
usflu.genind <- DNAbin2genind(x=usflu.dna, pop=usflu.annot[["year"]])
usflu.genind # Check it
# Read sequence data in NEXUS (similar to reading FASTA)
?read.nexus.data

# Importing DNA sequences from GeneBank - according to sequence ID, data from https://www.ncbi.nlm.nih.gov/popset/22854787
gunnera.dna <- read.GenBank(access.nb=c("AF447749.1", "AF447748.1", "AF447747.1", "AF447746.1", "AF447745.1", "AF447744.1", "AF447743.1", "AF447742.1", "AF447741.1", "AF447740.1", "AF447739.1", "AF447738.1", "AF447737.1", "AF447736.1", "AF447735.1", "AF447734.1", "AF447733.1", "AF447732.1", "AF447731.1", "AF447730.1", "AF447729.1", "AF447728.1"))
gunnera.dna
class(gunnera.dna)

# Query on-line sequence databases
library(seqinr)
choosebank() # Genetic banks available for seqinr
choosebank("embl", timeout=20) # Choose some bank
# Query selected database - there are much possibilities
?query # See how to construct the query
nothofagus <- query(listname="nothofagus", query="SP=Nothofagus AND K=rbcL", verbose=TRUE)
# See the sequences information
nothofagus$req
# Get the sequences as a list
nothofagus.sequences <- getSequence(nothofagus$req)
# See sequences
nothofagus.sequences
# Get annotations
nothofagus.annot <- getAnnot(nothofagus$req)
nothofagus.annot
# Close the bank when work is over
closebank()
# Convert sequences from a list to DNAbin (functions as.DNAbin*)
nothofagus.dna <- as.DNAbin.list(nothofagus.sequences)
nothofagus.dna # See it

# Query NCBI databases
library(rentrez)
entrez_dbs() # Genetic banks available for rentrez
entrez_db_summary("nucleotide") # Brief description of the selected database
entrez_db_searchable("nucleotide") # Set of search terms that can used with this database
# Search the database and get IDs of matched records
nothofagus.search <- entrez_search(db="nucleotide", term="Nothofagus[ORGN] AND rbcL[GENE]")
nothofagus.search
# Fetch desired records according to their IDs
nothofagus.fasta <- entrez_fetch(db="nucleotide", id=nothofagus.search[["ids"]], rettype="fasta")
nothofagus.fasta
# Conversion into DNAbin requires saving to disk as FASTA and then loading it
write(x=nothofagus.fasta, file="nothofagus.fasta")
nothofagus.dna <- read.dna(file="nothofagus.fasta", format="fasta")
nothofagus.dna

# Importing SNP
?read.PLINK # See it for available options
# Extract SNP from FASTA alignment
usflu.genlight <- fasta2genlight(file="https://adegenet.r-forge.r-project.org/files/usflu.fasta", quiet=FALSE, saveNbAlleles=TRUE)
# See result
usflu.genlight
# If it crashes (on Windows), try to add parameter "parallel=FALSE"
# Function has several options to speed up reading
?fasta2genlight

# Read SNP in Adegent format - check Adegenet tutorial first
?read.snp

# Checking SNPs
# Position of polymorphism within alignment - snpposi.plot() requires input data in form of matrix
snpposi.plot(x=as.matrix(usflu.dna), codon=FALSE)
# Position of polymorphism within alignment - differentiating codons
snpposi.plot(as.matrix(usflu.dna))
# When converting to genind object, only polymorphic loci are kept - threshold for polymorphism can be arbitrary (polyThres=...)
usflu.genind2 <- DNAbin2genind(x=usflu.dna, polyThres=0.01)
usflu.genind2 # See it
# Test is distribution of SNP is random (1000 permutations)
snpposi.test(x=as.matrix(usflu.dna))

# Check sequences
# Nucleotide diversity
pegas::nuc.div(x=usflu.dna)
# Base frequencies
ape::base.freq(x=usflu.dna)
# GC content
ape::GC.content(x=usflu.dna)
# Number of times any dimer/trimer/etc oligomers occur in a sequence
seqinr::count(seq=as.character.DNAbin(gunnera.dna[["AF447749.1"]]), wordsize=3)
# View sequences - all must be of the same length
# Function "image.DNAbin" requires as input matrix, so that sequences must be of same length (aligned)
image.DNAbin(x=usflu.dna)
# Sequences must be of same length - as.matrix.DNAbin() can help
image.DNAbin(x=as.matrix.DNAbin(usflu.dna))

# VCF
# Required library
library(vcfR)
# Read input file
# Download input file from https://soubory.trapa.cz/rcourse/arabidopsis.vcf.gz
arabidopsis.vcf <- read.vcfR(file=file.choose()) # Pick up downloaded file 'arabidopsis.vcf.gz' from the disc
# Or directly load remote file
arabidopsis.vcf <- read.vcfR(file="https://soubory.trapa.cz/rcourse/arabidopsis.vcf.gz")
# It returns object of class vcfR-class
?read.vcfR # See more import options
# Another option from package pegas returns list of objects loci and data.frame
?pegas::read.vcf

# Explore VCF
arabidopsis.vcf
head(arabidopsis.vcf)
arabidopsis.vcf@fix[1:10,1:5]
strwrap(arabidopsis.vcf@meta[1:21])
queryMETA(x=arabidopsis.vcf)
queryMETA(x=arabidopsis.vcf, element="DP")
queryMETA(x=arabidopsis.vcf, element="FORMAT.+DP")
queryMETA(x=arabidopsis.vcf, element="FORMAT=<ID=DP")
head(x=getFIX(x=arabidopsis.vcf))
head(x=is.polymorphic(x=arabidopsis.vcf, na.omit=TRUE))
head(x=is.biallelic(x=arabidopsis.vcf))
arabidopsis.vcf@gt[1:10, 1:4]
# See description of depth of coverage (DP) slot
strwrap(x=grep(pattern="ID=DP,", x=arabidopsis.vcf@meta, value=TRUE))

# Extract the depth of coverage (DP)
arabidopsis.vcf.dp <- extract.gt(x=arabidopsis.vcf, element="DP", as.numeric=TRUE) # GT:GQ:DP:HQ
# See it
dim(arabidopsis.vcf.dp)
head(arabidopsis.vcf.dp)
# Boxplot of DP
boxplot(x=arabidopsis.vcf.dp, col="#808080", ylab="Depth of coverage", las=3)
title("DP per specimen")
abline(h=seq(from=0, to=90, by=10), col="#b3b3b3")
# Bar plot of mean DP
barplot(apply(X=arabidopsis.vcf.dp, MARGIN=2, FUN=mean, na.rm=TRUE), las=3)
title("Mean DP per specimen")
abline(h=seq(from=0, to=60, by=10), col="#b3b3b3")
# Heatmap of DP (subset)
heatmap.bp(x=arabidopsis.vcf.dp[1:100,1:100], col.ramp=rainbow(n=100, start=0.1)) # Subset - only first 100 loci and individuals
title("DP per specimens and loci")

# Extract the genotype quality (GQ)
arabidopsis.vcf.gq <- extract.gt(x=arabidopsis.vcf, element="GQ", as.numeric=TRUE)
dim(arabidopsis.vcf.gq)
arabidopsis.vcf.gq[1:10,1:10]
# Heatmap of GQ (subset)
heatmap.bp(x=arabidopsis.vcf.gq[1:100,1:100])
# Bar plot of mean GQ
barplot(apply(X=arabidopsis.vcf.gq, MARGIN=2, FUN=mean, na.rm=TRUE), las=3)
abline(h=seq(from=0, to=90, by=5), col="grey")
# Boxplot of GQ
boxplot(arabidopsis.vcf.gq, las=2, main="Genotype Quality (GQ)")
abline(h=seq(from=0, to=100, by=5), col="grey")

# Missing data
# Extract information about missing data
arabidopsis.vcf.miss <- apply(X=arabidopsis.vcf.dp, MARGIN=2, FUN=function(x){sum(is.na(x))})
arabidopsis.vcf.miss <- arabidopsis.vcf.miss/nrow(arabidopsis.vcf)
# Bar plot of missing data
barplot(height=arabidopsis.vcf.miss, ylab="Percentage of missing data", las=2)
abline(h=seq(from=0, to=1, by=0.05), col="grey")
# Histogram of frequencies of missing data
arabidopsis.vcf.missg <- apply(X=arabidopsis.vcf.dp, MARGIN=1, FUN=function(x){sum(is.na(x))})
arabidopsis.vcf.missg <- arabidopsis.vcf.miss/ncol(arabidopsis.vcf@gt[,-1])
hist(x=arabidopsis.vcf.missg, xlab="Missingness (%)")
abline(h=seq(from=0, to=350, by=25), col="grey")

# Remove indels
arabidopsis.vcf <- extract.indels(x=arabidopsis.vcf)
# Remove non-biallelic loci
arabidopsis.vcf <- arabidopsis.vcf[is.biallelic(x=arabidopsis.vcf),]
# See result
arabidopsis.vcf

# ChromR - filtration of VCF
# Loading reference sequence
# Download https://soubory.trapa.cz/rcourse/alygenomes.fasta into working directory
arabidopsis.dna <- read.dna(file="alygenomes.fasta", format="fasta")
arabidopsis.dna
# Conversion into chromosome object (ChromR)
arabidopsis.chrom <- create.chromR(vcf=arabidopsis.vcf, name="Arabidopsis arenosa", seq=arabidopsis.dna)
arabidopsis.chrom
plot(arabidopsis.chrom)
# Masking sites with too low/high DP and/or MQ
arabidopsis.chrom.mask <- masker(x=arabidopsis.chrom, min_QUAL=1, min_DP=8, max_DP=5000, min_MQ=40, max_MQ=200)
arabidopsis.chrom.mask
plot(arabidopsis.chrom.mask)
variant.table(arabidopsis.chrom.mask)
# Saving mask into new object
arabidopsis.chrom.fin <- proc.chromR(x=arabidopsis.chrom.mask)
arabidopsis.chrom.fin
chromoqc(chrom=arabidopsis.chrom.fin) # The plot is bit empty as we have only single gene

# Convert VCF into various objects for later processing
# Genind
arabidopsis.genind <- vcfR2genind(x=arabidopsis.chrom.fin@vcf, ploidy=4)
# Check it
arabidopsis.genind
nInd(arabidopsis.genind)
indNames(arabidopsis.genind)
nLoc(arabidopsis.genind)
locNames(arabidopsis.genind)
# Genlight (suitable for huge data, not required now)
# Note that it introduces a lot of missing data due to variable ploidies
arabidopsis.genlight <- vcfR2genlight(x=arabidopsis.chrom.fin@vcf, n.cores=1) # On Linux/macOS and with large data use higher n.cores
# Check it
warnings() # See errors - due to missing data when handling tetraploids together with diploids
arabidopsis.genlight
# Loci
arabidopsis.loci <- vcfR2loci(x=arabidopsis.chrom.fin)
# Check it
arabidopsis.loci
print(x=arabidopsis.loci, details=TRUE)
# DNAbin
?vcfR2DNAbin # There are various options how to process variants in VCF
arabidopsis.dnabin <- vcfR2DNAbin(x=arabidopsis.chrom.fin, consensus=FALSE, extract.haps=TRUE, unphased_as_NA=FALSE, asterisk_as_del=FALSE)
# Check it
arabidopsis.dnabin
dim(arabidopsis.dnabin)
as.character.DNAbin(arabidopsis.dnabin[1:15,1:12])
image.DNAbin(arabidopsis.dnabin)
snpposi.plot.DNAbin(arabidopsis.dnabin)
snpposi.test.DNAbin(arabidopsis.dnabin)

## Exporting data
# Convert genind into DF using genind2df()
hauss.df <- genind2df(x=hauss.genind, pop=NULL, sep="/", usepop=TRUE, oneColPerAll=FALSE)
# Save microsatellites to disk - check settings of write.table
write.table(x=hauss.df, file="haussdata.txt", quote=FALSE, sep="\t", na="NA", dec=".", row.names=TRUE, col.names=TRUE)
# Export of DNA sequences into FASTA format
write.dna(x=usflu.dna, file="usflu.fasta", format="fasta", append=FALSE, nbcol=6)
seqinr::write.fasta(sequences=as.character.DNAbin(gunnera.dna), names=names(gunnera.dna), file.out="gunnera.fasta", open="w")
# Export DNA sequnces as NEXUS
write.nexus.data(x=gunnera.dna, file="gunnera.nexus", format="dna")
# Export VCF
vcfR::write.vcf(x=arabidopsis.vcf, file="arabidopsis.vcf.gz")
# Export trees (objects of class phylo)
write.tree(phy=hauss.nj.bruvo, file="haussknechtii.nwk") # Writes tree(s) in NEWICK format
write.nexus(hauss.nj.bruvo, file="haussknechtii.nexus") # Writes tree(s) in NEXUS format

## Multiple sequence alignment

# Libraries
library(ape)
library(ips)

# Requires path to MAFFT binary - set it according to your installation
# Read ?mafft and mafft's documentation
gunnera.mafft <- ips::mafft(x=gunnera.dna, method="localpair", maxiterate=100, options="--adjustdirection", exec="/usr/bin/mafft") # Change "exec" to fit your path to mafft (on Windows point to mafft.bat)!
gunnera.mafft # See results, compare with 'gunnera'
class(gunnera.mafft)
image.DNAbin(gunnera.mafft) # Plot the alignment

# Multiple sequence alignments using clustal, muscle and t-coffee are available in package ape
# read ?clustal and documentation of Clustal and Muscle to set correct parameters
gunnera.clustal <- ape::clustal(x=gunnera.dna, pw.gapopen=10, pw.gapext=0.1, gapopen=10, gapext=0.2, exec="/usr/bin/clustalw2", quiet=FALSE, original.ordering=TRUE) # Change "exec" to fit your path to clustal!
gunnera.clustal # See results
class(gunnera.clustal)
image.DNAbin(gunnera.clustal) # Plot the alignment
gunnera.muscle <- ape::muscle(x=gunnera.dna, exec="muscle", quiet=FALSE, original.ordering=TRUE) # Change "exec" to fit your path to muscle!
gunnera.muscle # See results
class(gunnera.muscle)
image.DNAbin(gunnera.muscle) # Plot the alignment
# See options in muscle package
?muscle::muscle

# Remove gaps from alignment - destroy it
gunnera.nogaps <- del.gaps(gunnera.muscle)
?del.gaps # See for details

# Align multiple genes
# Create a list of DNAbin objects to process
multialign <- list(gunnera.dna, usflu.dna, usflu.dna2)
# See it
multialign
class(multialign)
lapply(X=multialign, FUN=class)
# Do the alignment
multialign.aln <- lapply(X=multialign, FUN=ips::mafft, method="localpair", maxiterate=100, exec="/usr/bin/mafft") # Change "exec" to fit your path to mafft (on Windows point to mafft.bat)!
# See result
multialign.aln
multialign.aln[[1]]
lapply(X=multialign.aln, FUN=class)
# Do the same in parallel (mclapply do the tasks in parallel, not one-by-one like lapply)
library(parallel)
multialign.aln2 <- mclapply(X=multialign, FUN=ape::muscle, exec="muscle", quiet=FALSE, original.ordering=TRUE) # Change "path" to fit your path to muscle!
# mclapply() relies on forking and hence is not available on Windows unless "mc.cores=1"
# See result
multialign.aln2
lapply(X=multialign.aln2, FUN=class)
?mclapply # See more options
?clusterApply # See more options (parLapply should work on Windows)

# Plotting alignment
image.DNAbin(x=gunnera.mafft)
# Check the alignment
checkAlignment(x=usflu.dna, check.gaps=TRUE, plot=TRUE, what=1:4)
checkAlignment(x=as.matrix.DNAbin(x=gunnera.clustal), check.gaps=TRUE, plot=TRUE, what=1:4)
?checkAlignment # See details
# DNAbin can be technically list or matrix - some functions require list, some matrix, some can handle both - check manual and if needed, use
as.matrix.DNAbin()
as.list.DNAbin()
# Matrix makes sense only for alignments, list for any import (sequences do no have to have same lengths)

# Delete all columns/rows containing only gaps or missing data (N, ?, -)
gunnera.mafft <- deleteEmptyCells(DNAbin=gunnera.mafft)
?ips::deleteEmptyCells # See help page for details
?phyloch::delete.empty.cells # See help page for details
# Delete all columns containing at least 25% of gaps
gunnera.mafft.ng <- deleteGaps(x=gunnera.mafft, gap.max=nrow(gunnera.mafft)/4)
gunnera.mafft.ng
# Do not confuse with function delete.gaps() from phyloch package
# See of settings of "nmax" value - threshold for gap deletion
?deleteGaps # "nmax=0" deletes all columns with any gap
multialign.aln.ng <- lapply(X=multialign.aln, FUN=deleteGaps, gap.max=5)
multialign.aln.ng
# Delete every line (sample) containing at least 20% of missing data
gunnera.mafft.ng <- del.rowgapsonly(x=gunnera.mafft.ng, threshold=0.2, freq.only=FALSE)
gunnera.mafft.ng
?ape::del.rowgapsonly # See help page for details
# Delete every alignment position having at least 20% of missing data
gunnera.mafft.ng <- del.colgapsonly(x=gunnera.mafft.ng, threshold=0.2, freq.only=FALSE)
gunnera.mafft.ng
?ape::del.colgapsonly # See help page for details
# Display the result
image.DNAbin(x=gunnera.mafft.ng)
lapply(X=multialign.aln.ng, FUN=image.DNAbin)

## Descriptive statistics

# Load needed libraries
library(ape) # Analysis of phylogenetics and evolution
library(ade4) # Analysis of ecological data, multivariate methods
library(adegenet) # Exploratory analysis of genetic and genomic data
library(pegas) # Population and evolutionary genetics
library(poppr) # Population genetic analysis, including populations with mixed reproduction
library(hierfstat) # Hierarchical F-statistics
library(corrplot) # Visualization of correlation matrix
library(StAMPP) # Statistical analysis of mixed ploidy populations
library(philentropy) # Various genetic distances

# Get summary - names and sizes of populations, heterozygosity, some info about loci
hauss.summ <- summary(hauss.genind)
hauss.summ
# Plot expected vs. observed heterozygosity - it looks like big difference
plot(x=hauss.summ$Hexp, y=hauss.summ$Hobs, main="Observed vs expected heterozygosity", xlab="Expected heterozygosity", ylab="Observed heterozygosity")
abline(0, 1, col="red")
# Test normality of the vector - can we use bartlett.test() or t.test(), or do we have to use weaker kruskal.test() or wilcox.test()?
shapiro.test(hauss.summ$Hexp)
# Bartlett's K-squared test of difference between observed and expected heterozygosity - not significant
bartlett.test(list(hauss.summ$Hexp, hauss.summ$Hobs))
# T-test of difference between observed and expected heterozygosity - strongly significant
t.test(x=hauss.summ$Hexp, y=hauss.summ$Hobs, paired=TRUE, var.equal=TRUE)
# Create pane with some information
par(mfrow=c(2, 2)) # Divide graphical devices into 4 smaller spaces
# Plot alleles number vs. population sizes
plot(x=hauss.summ$n.by.pop, y=hauss.summ$pop.nall, xlab="Populations sample size", ylab="Number of alleles", main="Alleles numbers and sample sizes", col="red", pch=20)
# Add text description to the point
text(x=hauss.summ$n.by.pop, y=hauss.summ$pop.nall, lab=names(hauss.summ$n.by.pop), cex=1.5)
# Barplots of various data
barplot(height=hauss.summ$loc.n.all, ylab="Number of alleles", main="Number of alleles per locus", las=3)
barplot(height=hauss.summ$Hexp-hauss.summ$Hobs, main="Heterozygosity: expected-observed", ylab="Hexp - Hobs", las=3)
barplot(height=hauss.summ[["n.by.pop"]], main="Sample sizes per population", ylab="Number of genotypes", las=3)
dev.off() # Closes graphical device - otherwise following graphs would still be divided into 4 parts

# poppr returns various population statistics for populations
?poppr # See details
poppr(dat=hauss.genind, total=TRUE, sample=1000, method=4, missing="geno", cutoff=0.15, quiet=FALSE, clonecorrect=FALSE, plot=TRUE, index="rbarD", minsamp=1, legend=TRUE)
# Algorithms and equations utilized in poppr
vignette("algo", package="poppr")

# Departure from HWE
# According to loci
hauss.hwe.test <- hw.test(x=hauss.loci, B=1000)
# Results per-locus
haussknechtii.hwe.test
# Summary across all loci
summary(hauss.hwe.test)
# According to populations
# Separate genind object into list of genind objects for individual populations
hauss.pops <- seppop(hauss.genind)
hauss.pops
# Convert genind back to loci (list of loci objects according to populations)
hauss.pops.loci <- lapply(X=hauss.pops, FUN=genind2loci)
# Calculate the results per populations
lapply(X=hauss.pops.loci, FUN=hw.test, B=1000)

# FST
# Fit, Fst and Fis for each locus
# For Fst, fstat and theta.msat the loci object must contain population column
Fst(x=hauss.loci, pop=1) # Calculation per locus
summary(Fst(x=hauss.loci, pop=1)) # Summary across loci
# Nei's pairwise Fst between all pairs of populations.
pairwise.neifst(dat=genind2hierfstat(dat=hauss.genind))
heatmap(x=pairwise.neifst(dat=genind2hierfstat(dat=hauss.genind)), Rowv=NA, Colv=NA) # Simple visualization
# For mixed ploidy data sets
# stamppFst requires population factor in genlight (here, population code consists of first three letters of individual's name)
indNames(arabidopsis.genlight)
# Population code consists of first three letters of individual's name - extract the population name part
substr(x=indNames(arabidopsis.genlight), start=1, stop=3)
pop(arabidopsis.genlight) <- substr(x=indNames(arabidopsis.genlight), start=1, stop=3)
pop(arabidopsis.genlight) # Check it
popNames(arabidopsis.genlight)
arabidopsis.fst <- StAMPP::stamppFst(geno=arabidopsis.genlight, nboots=100, percent=95, nclusters=1)
# For large data use higher nclusters to parallelize calculations
?StAMPP::stamppFst # See details
# Matrix of Fst among populations
arabidopsis.fst[["Fsts"]]
# Matrix of P values
arabidopsis.fst[["Pvalues"]]
# Save results - open in spreadsheet (e.g. LibreOffice Calc)
write.table(x=arabidopsis.fst[["Fsts"]], file="arabidopsis_fst.tsv", quote=FALSE, sep="\t")
# Correlation plot of pairwise Fst
corrplot(corr=arabidopsis.fst[["Fsts"]], method="circle", type="lower", col=funky(15), title="Correlation matrix of Fst among populations", is.corr=FALSE, diag=FALSE, outline=TRUE, order="alphabet", tl.pos="lt", tl.col="black")
?corrplot # See for more options
# Display in similar way also another Fst tables

# TODO Estimates of population theta according to Kimmel et al. 1998
# theta.msat(hauss.loci)

## Multi locus genotypes
# Total number of MLGs (simple value)
mlg(gid=hauss.genind, quiet=FALSE)
# MLGs shared among populations (a list)
mlg.crosspop(gid=hauss.genind, df=TRUE, quiet=FALSE)
# Detailed view on distribution of MLGs into populations (table and/or plot)
mlg.table(gid=hauss.genind, bar=TRUE, total=TRUE, quiet=FALSE)
# Which individual belong to which vector (two ways of display)
mlg.vector(hauss.genind)
mlg.id(hauss.genind)

## Inbreeding
# Load training data
data(microbov)
# Separate populations of Salers cattle
microbov.pops <- seppop(microbov)[["Salers"]]
microbov.pops # See it
# Calculate the inbreeding
?inbreeding # Check for more settings
microbov.inbr <- inbreeding(x=microbov.pops, N=100)
# Prepare population means for plotting
microbov.bar <- sapply(X=microbov.inbr, FUN=mean)
# Plot it
hist(x=microbov.bar, col="firebrick", main="Average inbreeding in Salers cattle")

## MSN based on Bruvo's distance
# See details and options...
?bruvo.msn
# Get the MSN
# Note SSRs repeats 'rep(2, 12)' - change according to your data
# Here 12 SSRS loci of length 2 bp
bruvo.msn(gid=hauss.genind, replen=rep(x=2, times=12), loss=TRUE, palette=rainbow, vertex.label ="inds", gscale=TRUE, wscale=TRUE, showplot=TRUE)
# For another data types
?msn.poppr
# Interactive creation of MSN
?imsn

## Genetic distances

# Simple dissimilarity distance matrix
?dist.gene # Details about methods of this distance constructions
hauss.dist <- dist.gene(x=hauss.genind@tab, method="pairwise")
hauss.dist

# Nei's distance (not Euclidean) for populations (other methods are available, see ?dist.genpop)
hauss.dist.pop <- dist.genpop(x=hauss.genpop, method=1, diag=TRUE, upper=TRUE)
hauss.dist.pop
# Test if it is Euclidean
is.euclid(hauss.dist.pop, plot=TRUE, print=TRUE, tol=1e-10)
# Turns to be Euclidean
hauss.dist.pop <- cailliez(distmat=hauss.dist.pop, print=FALSE, tol=1e-07, cor.zero=TRUE)
# Test if it is Euclidean
is.euclid(hauss.dist.pop, plot=TRUE, print=TRUE, tol=1e-10)
# Show it
hauss.dist.pop

# Bruvo's distances weighting SSRs repeats - take care about replen parameter - requires repetition length for every SSRs locus
# E.g. if having 5 SSRs with repeat lengths 2, 2, 3, 3 and 2 bp use: bruvo.dist(pop=... replen=c(2, 2, 3, 3, 2)...)
hauss.dist.bruvo <- bruvo.dist(pop=hauss.genind, replen=rep(x=2, times=12), loss=TRUE)
# Test if it is Euclidean
is.euclid(hauss.dist.bruvo, plot=TRUE, print=TRUE, tol=1e-10)
# Turn to be Euclidean
hauss.dist.bruvo <- cailliez(distmat=hauss.dist.bruvo, print=FALSE, tol=1e-07, cor.zero=TRUE)
# Test if it is Euclidean
is.euclid(hauss.dist.bruvo, plot=TRUE, print=TRUE, tol=1e-10)
# Show it
hauss.dist.bruvo

# Nei's distance (not Euclidean) for individuals (other methods are available, see ?nei.dist from poppr package)
hauss.dist.nei <- poppr::nei.dist(x=hauss.genind, warning=TRUE)
# Test if it is Euclidean
is.euclid(distmat=hauss.dist.nei, plot=TRUE, print=TRUE, tol=1e-10)
# Show it
hauss.dist.nei

# Dissimilarity matrix returns a distance reflecting the number of allelic differences between two individuals
hauss.dist.diss <- diss.dist(x=hauss.genind, percent=FALSE, mat=TRUE)
# Test if it is Euclidean
is.euclid(distmat=as.dist(m=hauss.dist.diss), plot=TRUE, print=TRUE, tol=1e-10)
# Show it
hauss.dist.diss

# Import custom distance matrix - distances.txt must exist of the disk
MyDistance <- read.csv("distances.txt", header=TRUE, sep="\t", dec=".", row.names=1) # Or so...
MyDistance <- as.dist(MyDistance)
class(MyDistance)
is.euclid(MyDistance, plot=TRUE, print=TRUE, tol = 1e-10)
dim(MyDistance)
MyDistance
MyDistance <- cailliez(MyDistance, print=TRUE, tol = 1e-10, cor.zero=TRUE)
is.euclid(MyDistance, plot=TRUE, print=TRUE, tol = 1e-10)
MyDistance

# Compare different distance matrices
# List of functions to be parsed to respective dist.* function
distances <- c("Nei", "Rogers", "Edwards", "Reynolds", "Prevosti")
# Calculate the distance matrices
dists <- lapply(distances, function(x) {
	DISTFUN <- match.fun(paste(tolower(x), "dist", sep="."))
	DISTFUN(hauss.genind.cor)
	})
# Add names for the distance names
names(dists) <- distances
# Add Bruvo distance
dists[["Bruvo"]] <- hauss.dist.bruvo
dists
# Split graphical device into 2 lines, 3 panes each
par(mfrow=c(2, 3))
# Calculate NJ and plot all trees
x <- lapply(names(dists), function(x) {
	plot(njs(dists[[x]]), main=x, type="unrooted")
	add.scale.bar(lcol="red", length=0.1)
	})
dev.off() # Close graphical device to reset settings
# See details of distance methods in package poppr
vignette("algo", package="poppr")

# Selecting model for calculating distances among DNA sequences
library(phangorn) # Load needed library
?modelTest # Plenty of options, including possible parallelization
usflu.models <- modelTest(object=as.phyDat(usflu.dna)) # Test models
# Sort results according to AIC or BIC - the lower the better fit
head(usflu.models[order(usflu.models[["AIC"]]),], n=15L)
# Same for Gunnera dataset
gunnera.models <- modelTest(object=as.phyDat(gunnera.mafft.ng))
head(gunnera.models[order(gunnera.models[["AIC"]]),], n=15L)

# DNA distances - there are various models available
?dist.dna
usflu.dist <- dist.dna(x=usflu.dna, model="F84")
usflu.dist
class(usflu.dist)
dim(as.matrix(usflu.dist))
gunnera.dist <- dist.dna(x=gunnera.mafft.ng, model="TN93")
gunnera.dist
class(gunnera.dist)
dim(as.matrix(gunnera.dist))

# Distances and genlight object. Pairwise genetic distances for each data block (genlight objects with whole genome data) - sensitive to missing data
usflu.dists.l <- seploc(usflu.genlight, n.block=10, parallel=FALSE)
class(usflu.dists.l)
usflu.dists <- lapply(X=usflu.dists.l, FUN=function(DDD) dist(as.matrix(DDD)))
class(usflu.dists)
names(usflu.dists)
class(usflu.dists[[1]])
usflu.distr <- Reduce(f="+", x=usflu.dists)
class(usflu.distr)
usflu.distr
# It is possible to use just basic dist function on whole genlight object (might require a lot of RAM)
usflu.distg <- dist.gene(as.matrix(usflu.genlight))

# Distances in mixed-ploidy data sets
# stamppNeisD requires population factor in genlight
# Nei's 1972 distance between individuals (use pop=TRUE to calculate among populations)
arabidopsis.dist <- stamppNeisD(geno=arabidopsis.genlight, pop=FALSE)
# Check it
head(arabidopsis.dist)
dim(arabidopsis.dist)
class(arabidopsis.dist)
# The same on population level
arabidopsis.dist.pop <- stamppNeisD(geno=arabidopsis.genlight, pop=TRUE)
# Check it
head(arabidopsis.dist.pop)
dim(arabidopsis.dist.pop)
class(arabidopsis.dist.pop)
# Export the distance matrix as Phylip format for usage in external software (e.g. SplitsTree)
stamppPhylip(distance.mat=arabidopsis.dist, file="arabidopsis_dist.txt")
# Genomic relationship matrix
?stamppGmatrix # Method details
arabidopsis.genomat <- stamppGmatrix(geno=arabidopsis.genlight)
# Check it
head(arabidopsis.genomat)
dim(arabidopsis.genomat)
class(arabidopsis.genomat)

# Over 40 distances from philentropy package
getDistMethods() # See available distances
?distance # See details of distances
# Calculate e.g. Jaccard index for AFLP data
amara.jac <- distance(x=amara.aflp[,2:30], method="jaccard") # amara.aflp has 30 columns (see dim(amara.aflp), column 1 contains names (head(amara.aflp)))
# See result
class(amara.jac)
amara.jac
# Make it distance matrix
amara.jac <- as.dist(m=amara.jac, diag=TRUE, upper=TRUE)
amara.jac

# Visualize pairwise genetic similarities
# table.paint() requires data frame, dist can't be directly converted to DF
table.paint(df=as.data.frame(as.matrix(usflu.dist)), cleg=0, clabel.row=0.5, clabel.col=0.5)
# Same visualization, colored
# heatmap() reorders values
heatmap(x=as.matrix(usflu.dist), Rowv=NA, Colv=NA, symm=TRUE)
# Another possibility is to use function corrplot() from package corrplot

## Heatmaps
# Based on various distances
heatmap(as.matrix(hauss.dist), symm=TRUE, labRow=rownames(as.matrix(hauss.dist.bruvo)), labCol=colnames(as.matrix(hauss.dist.bruvo))) # hauss.dist doesn't contain names of individuals - we have to add them here
heatmap(as.matrix(hauss.dist.pop), symm=TRUE)
heatmap(as.matrix(hauss.dist.bruvo), symm=TRUE)
heatmap(as.matrix(hauss.dist.diss), symm=TRUE)
?heatmap # Other options

## AMOVA
# From package pegas (doesn't directly show percentage of variance)
hauss.pop <- pop(hauss.genind)
hauss.amova <- pegas::amova(hauss.dist~hauss.pop, data=NULL, nperm=1000, is.squared=FALSE)
# See results
hauss.amova
# For more complicated hierarchy
?poppr::poppr.amova
# For mixed-ploidy dat sets
?StAMPP::stamppAmova

## Hierarchical clustering
# According to distance used
?hclust # How to use hierarchical clustering
plot(hclust(d=hauss.dist, method="complete"))
plot(hclust(d=hauss.dist.pop, method="complete"))
plot(hclust(d=hauss.dist.bruvo, method="complete"))

## UPGMA
# Calculate it
# Saving as phylo object (and not hclust) gives more possibilities for further plotting and manipulations
usflu.upgma <- as.phylo(hclust(d=usflu.dist, method="average"))
plot.phylo(x=usflu.upgma, cex=0.75)
title("UPGMA tree")
# Test quality - tests correlation of original distance in the matrix and reconstructed distance from hclust object
plot(x=as.vector(usflu.dist), y=as.vector(as.dist(cophenetic(usflu.upgma))), xlab="Original pairwise distances", ylab="Pairwise distances on the tree", main="Is UPGMA appropriate?", pch=20, col=transp(col="black", alpha=0.1), cex=2)
abline(lm(as.vector(as.dist(cophenetic(usflu.upgma)))~as.vector(usflu.dist)), col="red")

## NJ

# Calculates the tree (try with various distances)
hauss.nj <- nj(hauss.dist)

# Test tree quality - plot original vs. reconstructed distance
plot(as.vector(hauss.dist), as.vector(as.dist(cophenetic(hauss.nj))), xlab="Original distance", ylab="Reconstructed distance")
abline(lm(as.vector(hauss.dist) ~ as.vector(as.dist(cophenetic(hauss.nj)))), col="red")
cor.test(x=as.vector(hauss.dist), y=as.vector(as.dist(cophenetic(hauss.nj))), alternative="two.sided") # Testing the correlation
summary(lm(as.vector(hauss.dist) ~ as.vector(as.dist(cophenetic(hauss.nj))))) # Prints summary text

# Plot a basic tree - see ?plot.phylo for details
plot.phylo(x=hauss.nj, type="phylogram")
plot.phylo(x=hauss.nj, type="cladogram", edge.width=2)
plot.phylo(x=hauss.nj, type="fan", edge.width=2, edge.lty=2)
plot.phylo(x=hauss.nj, type="radial", edge.color="red", edge.width=2, edge.lty=3, cex=2)

# Bootstrap

# boot.phylo() resamples all columns - remove population column first
hauss.loci.nopop <- hauss.loci
hauss.loci.nopop[["population"]] <- NULL
hauss.boot <- boot.phylo(phy=hauss.nj, x=hauss.loci.nopop, FUN=function(XXX) nj(dist.gene(loci2genind(XXX)@tab, method="pairwise")), B=1000)
# boot.phylo returns NUMBER of replicates - NO PERCENTAGE
# Plot the tree
plot.phylo(x=hauss.nj, type="unrooted", main="Neighbour-Joining tree")
# Labels for nodes - bootstrap - see ?nodelabels for graphical settings
nodelabels(text=round(hauss.boot/10))

# Another bootstrap possibility
hauss.aboot <- aboot(x=hauss.genind, tree=nj, distance=nei.dist, sample=100) # Bootstrap values are in slot node.label
# Plot the tree, explicitly display node labels
plot.phylo(x=hauss.aboot, show.node.label=TRUE)
?aboot # See details...
?plot.phylo

# Plot a nice tree with colored tips
plot.phylo(x=hauss.nj, type="unrooted", show.tip=FALSE, edge.width=3, main="Neighbour-Joining tree")
# Labels for nodes - bootstrap - see ?nodelabels for graphical settings
nodelabels(text=round(hauss.boot/10))
# Coloured labels - creates vector of colors according to population information in genind object
nj.rainbow <- colorRampPalette(rainbow(length(popNames(hauss.genind))))
# Colored tips
tiplabels(text=indNames(hauss.genind), bg=fac2col(x=pop(hauss.genind), col.pal=nj.rainbow))

# Plot BW tree with tip symbols and legend
plot.phylo(x=hauss.nj, type="cladogram", show.tip=FALSE, edge.width=3, main="Neighbour-Joining tree")
# Add axis with distances
axisPhylo()
# From node labels let's remove unneeded frame
nodelabels(text=round(hauss.boot/10), frame="none", bg="white")
# As tip label we use only symbols - see ?points for graphical details
tiplabels(frame="none", pch=rep(x=0:4, times=c(13, 17, 2, 6, 9)), lwd=2, cex=2)
# Plot a legend explainindg symbols
legend(x="topleft", legend=c("He", "Oh", "Pr", "Ne", "Sk"), border="black", pch=0:4, pt.lwd=2, pt.cex=2, bty="o", bg="lightgrey", box.lwd=2, cex=1.2, title="Populations")
?plot.phylo
?nodelabels
?legend
?axisPhylo

# Bruvo's distance - NJ
hauss.nj.bruvo <- bruvo.boot(pop=hauss.genind, replen=rep(2, 12), sample=1000, tree="nj", showtree=TRUE, cutoff=1, quiet=FALSE)
plot.phylo(x=hauss.nj.bruvo, type="unrooted", show.tip=FALSE, edge.width=3, main="Neighbor-Joining tree")
# bruvo.boot() writes all needed information into resulting object, so there is no need for external bootstrap function
nodelabels(hauss.nj.bruvo[["node.labels"]]) # Note you can call node labels from phylo object as phylo$node.labels or phylo[["node.labels"]]
tiplabels(hauss.nj.bruvo[["tip.label"]], bg=fac2col(x=hauss.genind$pop, col.pal=nj.rainbow))

# Bruvo's distance - UPGMA
hauss.upgma <- bruvo.boot(pop=hauss.genind, replen=rep(2, 12), sample=1000, tree="upgma", showtree=TRUE, cutoff=1, quiet=FALSE)
plot.phylo(hauss.upgma, type="unrooted", show.tip=FALSE, edge.width=3, main="UPGMA tree")
nodelabels(hauss.upgma[["node.labels"]])
tiplabels(hauss.upgma[["tip.label"]], bg=fac2col(x=hauss.genind@pop, col.pal=nj.rainbow))

# Populations
# aboot() can use distances implemented in poppr:
?poppr::aboot
?poppr::nei.dist
# Calculations
hauss.nj.pop <- aboot(x=hauss.genpop, tree="nj", distance="nei.dist", sample=1000, showtree=FALSE)
# Information about result
print.phylo(hauss.nj.pop)
# Plot a tree
plot.phylo(x=hauss.nj.pop, type="radial", show.node.label=TRUE, cex=1.2, edge.width=3, main="Neighbor-Joining tree of populations")

# Make a tree based on DNA sequences
usflu.tree <- nj(X=usflu.dist)
# Plot it
plot.phylo(x=usflu.tree, type="unrooted", show.tip=FALSE)
title("Unrooted NJ tree")
# Coloured tips
usflu.pal <- colorRampPalette(topo.colors(length(levels(as.factor(usflu.annot[["year"]])))))
# Tip labels
tiplabels(text=usflu.annot$year, bg=num2col(usflu.annot$year, col.pal=usflu.pal), cex=0.75)
# Legend - describing years - pretty() automatically shows best values from given range, num2col() selects colors from color scale
legend(x="bottomright", fill=num2col(x=pretty(x=1993:2008, n=8), col.pal=usflu.pal), leg=pretty(x=1993:2008, n=8), ncol=1)

# Root the tree - "outgroup" is name of accession (in quotation marks) or number (position within phy object)
usflu.tree.rooted <- root.phylo(phy=usflu.tree, outgroup=1)
# Plot it
plot.phylo(x=usflu.tree.rooted, show.tip=FALSE, edge.width=2)
title("Rooted NJ tree")
# Labeling of tips
tiplabels(text=usflu.annot$year, bg=transp(num2col(x=usflu.annot$year, col.pal=usflu.pal), alpha=0.7), cex=0.75, fg="transparent")
# Add axis with phylogenetic distance
axisPhylo()
# Legend - describing years - pretty() automatically shows best values from given range, num2col() selects colors from color scale
legend(x="topright", fill=num2col(x=pretty(x=1993:2008, n=8), col.pal=usflu.pal), leg=pretty(x=1993:2008, n=8), ncol=1)

# Bootstrap rooted tree
# Calculate it
usflu.boot <- boot.phylo(phy=usflu.tree.rooted, x=usflu.dna, FUN=function(EEE) root.phylo(nj(dist.dna(EEE, model="TN93")), outgroup=1), B=1000)
# Plot the tree
plot.phylo(x=usflu.tree.rooted, show.tip=FALSE, edge.width=2)
title("NJ tree + bootstrap values")
tiplabels(frame="none", pch=20, col=transp(num2col(x=usflu.annot$year, col.pal=usflu.pal), alpha=0.7), cex=3.5, fg="transparent")
axisPhylo()
# Legend - describing years - pretty() automatically shows best values from given range, num2col() selects colors from color scale
legend(x="topright", fill=num2col(x=pretty(x=1993:2008, n=8), col.pal=usflu.pal), leg=pretty(x=1993:2008, n=8), ncol=1)
# Plots bootstrap support - note usflu.boot contains raw numbers - transform it into percent
nodelabels(text=round(usflu.boot/10), cex=0.75)

# Collapse branches with low bootstrap support
usflu.tree.usflu.na.density <- usflu.tree.rooted
# Determine branches with low support - note BS values are in raw numbers - use desired percentage with respect to number of bootstraps
usflu.tocollapse <- match(x=which(usflu.boot < 700)+length(usflu.tree.rooted$tip.label), table=usflu.tree.usflu.na.density$edge[,2])
# Set length of bad branches to zero
usflu.tree.usflu.na.density$edge.length[usflu.tocollapse] <- 0
# Create new tree
usflu.tree.collapsed <- di2multi(phy=usflu.tree.usflu.na.density, tol=0.00001)
# Plot the consensus tree
plot.phylo(x=usflu.tree.collapsed, show.tip=FALSE, edge.width=2)
title("NJ tree after collapsing weak nodes")
tiplabels(text=usflu.annot$year, bg=transp(num2col(x=usflu.annot$year, col.pal=usflu.pal), alpha=0.7), cex=0.5, fg="transparent")
axisPhylo()
legend(x="topright", fill=num2col(x=pretty(x=1993:2008, n=8), col.pal=usflu.pal), leg=pretty(x=1993:2008, n=8), ncol=1)

## Replacements of NJ
# Usage is same as nj() Check
?njs
?bionj
?bionjs
?NJ
?UNJ
?fastme

## PCoA

# Calculate it - based on various distance matrix
hauss.pcoa <- dudi.pco(d=dist.gene(x=scaleGen(x=hauss.genind, center=TRUE, scale=FALSE, truenames=TRUE), method="pairwise"), scannf=FALSE, nf=3)
# Basic display
s.label(dfxy=hauss.pcoa$li, clabel=0.75)
# To plot different axes use for example "...dfxy=hauss.pcoa$li[c(2, 3)]..."
# Add kernel density
s.kde2d(dfxy=hauss.pcoa$li, cpoint=0, add.plot=TRUE)
# Add histogram of Eigenvalues
add.scatter.eig(w=hauss.pcoa$eig, nf=3, xax=1, yax=2, posi="topleft", sub="Eigenvalues")
# Percentage of variance explained by each PC axis
100*hauss.pcoa$eig/sum(hauss.pcoa$eig)

# Colored display according to populations
hauss.pcoa.col <- rainbow(length(popNames(hauss.genind))) # Creates vector of colors according to populations
s.class(dfxy=hauss.pcoa$li, fac=pop(hauss.genind), col=hauss.pcoa.col)
add.scatter.eig(w=hauss.pcoa$eig, nf=3, xax=1, yax=2, posi="bottomleft", sub="Eigenvalues")
title("Principal Coordinates Analysis") # Adds title to the graph

# Based on Bruvo's distance
hauss.pcoa.bruvo <- dudi.pco(d=bruvo.dist(pop=hauss.genind, replen=rep(2, 12)), scannf=FALSE, nf=3)
s.class(dfxy=hauss.pcoa.bruvo$li, fac=pop(hauss.genind), col=hauss.pcoa.col)
add.scatter.eig(hauss.pcoa.bruvo$eig, posi="bottomright", 3, 1, 2)
# Another possibility for colored plot (see ?colorplot for details)
colorplot(xy=hauss.pcoa$li[c(1, 2)], X=hauss.pcoa$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title(main="PCoA, axes 1 and 2")
abline(v=0, h=0, col="grey", lty=2)

## SNP

# Plot of missing data (white) and number of 2nd alleles
glPlot(x=usflu.genlight, legend=TRUE, posi="topleft")
# Sum of the number of second allele in each SNP
usflu.freq <- glSum(usflu.genlight)
# Plot distribution of (second) allele frequencies
hist(x=usflu.freq, proba=TRUE, col="gold", xlab="Allele frequencies", main="Distribution of (second) allele frequencies")
lines(x=density(usflu.freq)$x, y=density(usflu.freq)$y*1.5, col="red", lwd=3 )
# Mean number of second allele in each SNP
usflu.mean <- glMean(usflu.genlight)
usflu.mean <- c(usflu.mean, 1-usflu.mean)
# Plot distribution of allele frequencies - play with parameters to get nice image
hist(x=usflu.mean, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies", main="Distribution of allele frequencies", nclass=20)
lines(x=density(usflu.mean, bw=0.05)$x, y=density(usflu.mean, bw=0.05)$y*2, lwd=3)

# Number of missing values in each locus
usflu.na.density <- density(glNA(usflu.genlight), bw=10) # Play with bw parameter to get optimal image
# Set range of xlim parameter from 0 to the length of original alignment
plot(x=usflu.na.density, type="n", xlab="Position in the alignment", main="Location of the missing values (NA)", xlim=c(0, 1701))
polygon(c(usflu.na.density$x, rev(usflu.na.density$x)), c(usflu.na.density$y, rep(0, length(usflu.na.density$x))), col=transp("blue", alpha=0.3))
points(glNA(usflu.genlight), rep(0, nLoc(usflu.genlight)), pch="|", cex=2, col="blue")

# PCA - select number of retained PC axes, about 10 here
usflu.pca <- glPca(x=usflu.genlight, center=TRUE, scale=FALSE, loadings=TRUE)
# Plot PCA
scatter.glPca(x=usflu.pca, posi="bottomright")
title("PCA of the US influenza data")
# Loading plot - contribution of variables to the pattern observed
loadingplot.glPca(x=usflu.pca)
# Coloured plot
colorplot(usflu.pca$scores, usflu.pca$scores, transp=TRUE, cex=4)
title("PCA of the US influenza data")
abline(h=0, v=0, col="grey")
add.scatter.eig(usflu.pca$eig[1:40], 2, 1, 2, posi="topright", inset=0.05, ratio=0.3)
# Calculate phylogenetic tree
usflu.tree.genlight <- nj(dist.gene(as.matrix(usflu.genlight)))
# Plot colored phylogenetic tree
plot.phylo(x=usflu.tree.genlight, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=num2col(usflu.annot[["year"]], col.pal=usflu.pal), cex=4)
title("NJ tree of the US influenza data")
# Add legend
legend(x="topright", fill=num2col(x=pretty(x=1993:2008, n=8), col.pal=usflu.pal), leg=pretty(x=1993:2008, n=8), ncol=1)

## Clustering analysis

# Needed library
library(adegenet)
# Tutorial, more information about DAPC
adegenetTutorial("dapc")
# Graphical web interface for DAPC
adegenetServer("DAPC")

## K-find
# Bayesian K-means clustering
# Retain all informative PC (here about 35). According to second graph select best K (here 2 or 3).
hauss.kfind <- find.clusters(x=hauss.genind, stat="BIC", choose.n.clust=TRUE, max.n.clust=10, n.iter=100000, n.start=100, scale=FALSE, truenames=TRUE)
# See results as text
table(pop(hauss.genind), hauss.kfind$grp)
hauss.kfind
# Graph showing table of original and inferred populations and assignment of individuals
table.value(df=table(pop(hauss.genind), hauss.kfind$grp), col.lab=paste("Inferred\ncluster", 1:length(hauss.kfind$size)), grid=TRUE)
# For K=3 - note parameters n.pca and n.clust - we just rerun the analysis and when results are stable, no problem here
hauss.kfind3 <- find.clusters(x=hauss.genind, n.pca=35, n.clust=3, stat="BIC", choose.n.clust=FALSE, n.iter=100000, n.start=100, scale=FALSE, truenames=TRUE)
# See results as text
table(pop(hauss.genind), hauss.kfind3$grp)
hauss.kfind3
# Graph showing table of original and inferred populations and assignment of individuals
table.value(df=table(pop(hauss.genind), hauss.kfind3$grp), col.lab=paste("Inferred\ncluster", 1:length(hauss.kfind3$size)), grid=TRUE)

## DAPC
# K=2
# Create DAPC
# Number of informative PC (Here 15, adegenet recommends < N/3). Select number of informative DA (here 1 - only one is available - it won't produce PCA graph).
hauss.dapc <- dapc(x=hauss.genind, pop=hauss.kfind$grp, center=TRUE, scale=FALSE, var.contrib=TRUE, pca.info=TRUE, truenames=TRUE)
# Information
hauss.dapc
# Density functions
scatter(x=hauss.dapc, xax=1, yax=1, main="DAPC", bg="white", solid=0.5, leg=TRUE, txt.leg=c("Group 1", "Group 2"), posi.leg="topright")
# Assignment of individuals to clusters
assignplot(x=hauss.dapc)
# Structure-like plot
compoplot(x=hauss.dapc, xlab="Individuals", leg=FALSE)
# Loadingplot - alleles the most adding to separation of individuals
loadingplot(x=hauss.dapc$var.contr)
# alfa-score - according to number of PC axis
optim.a.score(x=hauss.dapc)

# K=3
# Create DAPC
# Number of informative PC (Here 15, adegenet recommends < N/3). Select number of informative DA (here 2).
hauss.dapc3 <- dapc(x=hauss.genind, pop=hauss.kfind3$grp, center=TRUE, scale=FALSE, var.contrib=TRUE, pca.info=TRUE, truenames=TRUE)
# Information
hauss.dapc
# A la PCA graph
scatter(x=hauss.dapc3, main="DAPC, Taraxacum haussknechtii", bg="white", cex=3, clab=0, col=rainbow(3), posi.da="bottomleft", scree.pca=TRUE, posi.pca="bottomright", leg=TRUE, txt.leg=c("Group 1", "Group 2", "Group 3"), posi.leg="topleft")
# Same in BW
scatter(x=hauss.dapc3, main="DAPC, Taraxacum haussknechtii", bg="white", pch=c(15:17), cell=0, cstar=0, solid=1, cex=2.5, clab=0, col=gray.colors(3, start=0, end=0.8, gamma=2, alpha=0), posi.da="bottomleft", scree.pca=TRUE, posi.pca="bottomright", leg=TRUE, txt.leg=c("Group 1", "Group 2", "Group 3"), posi.leg="topleft")
# Density functions - only for first axis here!
scatter(x=hauss.dapc3, xax=1, yax=1, main="DAPC", bg="white", solid=0.5, leg=TRUE, txt.leg=c("Group 1", "Group 2", "Group 3"), posi.leg="topleft")
# Assignment of individuals to clusters
assignplot(hauss.dapc3)
# Structure-like plot
compoplot(hauss.dapc3, xlab="Individuals", leg=FALSE)
# Loadingplot - alleles the most adding to separation of individuals
loadingplot(hauss.dapc3$var.contr)
# alfa-score - according to number of PC axis
optim.a.score(hauss.dapc3)

## Moran's I
# Load required libraries
library(adegenet)
library(spdep)
# Creates connection network
hauss.connectivity <- chooseCN(xy=hauss.genind$other$xy, type=5, d1=0, d2=1, plot.nb=TRUE, result.type="listw", edit.nb=FALSE)
hauss.connectivity
# Test of Moran's I for 1st PCoA axis
moran.test(x=hauss.pcoa$li[,1], listw=hauss.connectivity, alternative="greater", randomisation=TRUE)
hauss.pcoa1.mctest <- moran.mc(x=hauss.pcoa[["li"]][,1], listw=hauss.connectivity, alternative="greater", nsim=1000)
hauss.pcoa1.mctest
# Plot the results
plot(hauss.pcoa1.mctest) # Plot of densitiy of permutations
moran.plot(x=hauss.pcoa$li[,1], listw=hauss.connectivity) # PC plot
# Test of Moran's I for 2nd PCoA axis
moran.test(x=hauss.pcoa$li[,2], listw=hauss.connectivity, alternative="greater", randomisation=TRUE)
hauss.pcoa2.mctest <- moran.mc(x=hauss.pcoa$li[,2], listw=hauss.connectivity, alternative="greater", nsim=1000)
hauss.pcoa2.mctest
plot(hauss.pcoa2.mctest) # Plot of densitiy of permutations
moran.plot(x=hauss.pcoa[["li"]][,2], listw=hauss.connectivity) # PC plot

## sPCA

# Explore various networks
data(rupica)
# Part of sPCA calculations are in adespatial package, load it
library(adespatial)
# Try various settings for chooseCN (type=X) - type 1-4 as there are identical coordinates (multiple sampling from same locality)
?chooseCN # See for more details - select the best "type" for your data
chooseCN(xy=rupica$other$xy, ask=TRUE, type=5/6/7, plot.nb=TRUE, edit.nb=FALSE, ...) # Play with settings little bit...

# Calculates sPCA - here is only 1 positive and 3 negative factors
hauss.spca <- spca(obj=hauss.genind, cn=hauss.connectivity, scale=TRUE, scannf=TRUE)
# Plot eigenvalues of sPCA - gobal vs. local structure
barplot(height=hauss.spca$eig, main="Eigenvalues of sPCA", col=spectral(length(hauss.spca[["eig"]])))
legend("topright", fill=spectral(2), leg=c("Global structures", "Local structures")) # Add legend
abline(h=0, col="grey") # Add line showing zero
# Information about sPCA
print.spca(hauss.spca)
# Summary of sPCA results
summary.spca(hauss.spca)
# Shows connectivity network, 3 different scores, barplot of eigenvalues and eigenvalues decomposition
plot.spca(hauss.spca)
# Display of scores in color canals (two analogous variants)
colorplot.spca(hauss.spca, cex=3)
title("sPCA - colorplot of PC 1 and 2 (lagged scores)", line=1, cex=1.5)
colorplot(x=hauss.genind$other$xy, hauss.spca$ls, axes=1:ncol(hauss.spca$li), cex=3)
title("sPCA - colorplot of PC 1 and 2 (lagged scores)", line=1, cex=1.5)
# Spatial and variance components of the eigenvalues
screeplot.spca(x=hauss.spca, main=NULL)

# Test if global/local structure is significant
hauss.spca.glo <- global.rtest(X=hauss.genind$tab, listw=hauss.spca$lw, nperm=999)
hauss.spca.glo
plot(hauss.spca.glo)
hauss.spca.loc <- local.rtest(X=hauss.genind$tab, listw=hauss.spca$lw, nperm=999)
hauss.spca.loc
plot(hauss.spca.loc)

# Map of genetic clines
library(akima) # This library is needed for some manipulation with coordinates
hauss.spca.temp <- interp(other(hauss.genind)$xy[,1], other(hauss.genind)$xy[,2], hauss.spca$ls[,1], xo=seq(min(other(hauss.genind)$xy[,1]), max(other(hauss.genind)$xy[,1]), le=200), yo=seq(min(other(hauss.genind)$xy[,2]), max(other(hauss.genind)$xy[,2]), le=200), duplicate="median")
# For 1st axis
image(x=hauss.spca.temp, col=spectral(100))
s.value(dfxy=hauss.genind$other$xy, z=hauss.pcoa$li[,1], add.p=TRUE, csize=0.5, sub="PCoA - first PC", csub=2, possub="topleft")
# For 2nd axis
image(x=hauss.spca.temp, col=spectral(100))
s.value(dfxy=hauss.genind$other$xy, z=hauss.pcoa[["li"]][,2], add.p=TRUE, csize=0.5, sub="PCoA - second PC", csub=2, possub="topleft")
# Interpolated lagged score on a map
hauss.spca.annot <- function() {
	title("sPCA - interpolated map of individual scores")
	points(other(hauss.genind)$xy[,1], other(hauss.genind)$xy[,2])
	}
filled.contour(hauss.spca.temp, color.pal=colorRampPalette(lightseasun(100)), pch=20, nlevels=100, key.title=title("Lagged\nscore 1"), plot.title=hauss.spca.annot())

# Loading plots - which alleles contribute the most?
hauss.spca.loadings <- hauss.spca$c1[,1]^2
names(hauss.spca.loadings) <- rownames(hauss.spca$c1)
loadingplot(x=hauss.spca.loadings, xlab="Alleles", ylab="Weight of the alleles", main="Contribution of alleles to the first sPCA axis")
boxplot(formula=hauss.spca.loadings~hauss.genind$loc.fac, las=3, ylab="Contribution", xlab="Marker", main="Contribution by markers into the first global score", col="gray")

## Mantel test - isolation by distance
# Geodesic distance
library(pegas)
# Geographical distance
hauss.gdist <- as.dist(m=geod(lon=hauss.genind$other$xy$lon, lat=hauss.genind$other$xy$lat), diag=TRUE, upper=TRUE)
# Mantel test
hauss.mantel <- mantel.randtest(m1=hauss.dist, m2=hauss.gdist, nrepet=1000)
hauss.mantel # See text output
plot(hauss.mantel, nclass=30)
# Libraries required by mantel.correlog:
library(permute)
library(lattice)
library(vegan)
# Different implementation of Mantel test testing distance classes
hauss.mantel.cor <- mantel.correlog(D.eco=hauss.dist, D.geo=hauss.gdist, XY=NULL, n.class=0, break.pts=NULL, cutoff=FALSE, r.type="pearson", nperm=1000, mult="holm", progressive=TRUE)
# See results for respective classes
hauss.mantel.cor
summary(hauss.mantel.cor)
# Plot it
plot(hauss.mantel.cor)

## Monmonier - genetic boundaries
# It requires every point to have unique coordinates (one can use jitter() or difference in scale of meters). Example here is on population level, which is not ideal.
# Calculates Monmonier's function (for threshold use 'd')
hauss.monmonier <- monmonier(xy=hauss.genpop$other$xy, dist=rogers.dist(x=hauss.genpop), cn=chooseCN(hauss.genpop$other$xy, ask=FALSE, type=6, k=2, plot.nb=FALSE, edit.nb=FALSE), nrun=1)
coords.monmonier(hauss.monmonier) # See result as text
# Plot genetic boundaries
plot.monmonier(hauss.monmonier, add.arrows=FALSE, bwd=10, sub="Monmonier plot", csub=2)
points(hauss.genpop$other$xy, cex=2.5, pch=20, col="red")
text(x=hauss.genpop$other$xy$lon, y=hauss.genpop$other$xy$lat, labels=popNames(hauss.genpop), cex=3)
legend("bottomright", leg="Genetic boundaries\namong populations")
# For plotting see
?points
?text
?legend

## Geneland
# Haploid and diploid codominant markers (microsattelites or SNPs)
# Load needed libraries
library(PBSmapping) # Required to transform coordinates
library(Geneland)
# Graphical interface available for Geneland - we will use only command line
Geneland.GUI()
# Loading and conversions of coordinates - Geneland requires specific coordinate space
hauss.geneland.coord <- as.matrix(hauss.coord) # hauss.cord is DF, we need just plain matrix
colnames(hauss.geneland.coord) <- c("X", "Y")
attr(hauss.geneland.coord, "projection") <- "LL"
attr(hauss.geneland.coord, "zone") <- NA
hauss.geneland.coord.utm <- convUL(hauss.geneland.coord)
dim(hauss.geneland.coord)
hauss.geneland.coord
dim(hauss.geneland.coord.utm)
hauss.geneland.coord.utm # Final coordinates
# Load data (only haploid or diploid data are supported) - only plain table with alleles
hauss.geneland.data <- read.table(file="https://soubory.trapa.cz/rcourse/haussknechtii_geneland.txt", na.string="-999", header=FALSE, sep="\t")
dim(hauss.geneland.data)
hauss.geneland.data
# Set number of independent runs
hauss.geneland.nrun <- 5
# Set length of burnin chain
hauss.geneland.burnin <- 100
# Set maximal K (number of populations)
hauss.geneland.maxpop <- 10
# In practice set much higher number of iterations (nit, millions), appropriate sampling (thinning, thousands) and longer burnin
# Number of iterations (MCMC steps)
hauss.geneland.nit <- 10000
# Sampling (thinning)
hauss.geneland.thinning <- 10
# FOR loop will run several independent runs and produce output maps of genetic clusters - outputs are written into subdirectory within geneland directory
for (hauss.geneland.irun in 1:hauss.geneland.nrun) {
	hauss.geneland.path.mcmc <- paste("geneland/", hauss.geneland.irun, "/", sep="") # paste is good especially for joining several text chains into one
	# On Windows, remove following line and create subdirectories from 1 to max K manually (creating subdirs in Windows in R is complicated)
	system(paste("mkdir ", hauss.geneland.path.mcmc)) # Creates subdirectory
	# Inferrence - MCMC chain - see ?MCMC for details
	MCMC(coordinates=hauss.geneland.coord.utm, geno.dip.codom=hauss.geneland.data, path.mcmc=hauss.geneland.path.mcmc, delta.coord=0.001, varnpop=TRUE, npopmin=1, npopmax=hauss.geneland.maxpop, nit=hauss.geneland.nit, thinning=hauss.geneland.thinning, freq.model="Uncorrelated", spatial=TRUE)
	# Post-process chains
	PostProcessChain(coordinates=hauss.geneland.coord.utm, path.mcmc=hauss.geneland.path.mcmc, nxdom=500, nydom=500, burnin=hauss.geneland.burnin)
	# Output
	# Simulated number of populations
	Plotnpop(path.mcmc=hauss.geneland.path.mcmc, printit=TRUE, file=paste(hauss.geneland.path.mcmc, "/geneland-number_of_clusters.pdf", sep=""), format="pdf", burnin=hauss.geneland.burnin)
	dev.off() # We must close graphical device manually
	# Map of estimated population membership
	PosteriorMode(coordinates=hauss.geneland.coord.utm, path.mcmc=hauss.geneland.path.mcmc, printit=TRUE, format="pdf", file=paste(hauss.geneland.path.mcmc,"/geneland-map.pdf", sep=""))
	dev.off() # We must close graphical device manually
	}
# Prepare list to record values of Fst for all runs
hauss.geneland.fstat <- list()
# Estimate Fst
for (hauss.geneland.irun in 1:hauss.geneland.nrun) {
	hauss.geneland.path.mcmc <- paste("geneland/", hauss.geneland.irun, "/", sep="")
	# F-statistics - Fis and Fst
	hauss.geneland.fstat[[hauss.geneland.irun]] <- Fstat.output(coordinates=hauss.geneland.coord.utm, genotypes=hauss.geneland.data, burnin=hauss.geneland.burnin, ploidy=2, path.mcmc=hauss.geneland.path.mcmc)
	}
# Print Fst output
hauss.geneland.fstat
# MCMC inference under the admixture model
for (hauss.geneland.irun in 1:hauss.geneland.nrun) {
	hauss.geneland.path.mcmc <- paste("geneland/", hauss.geneland.irun, "/", sep="")
	hauss.geneland.path.mcmc.adm <- paste(hauss.geneland.path.mcmc, "admixture", "/", sep="")
	# On Windows, remove following line of code and create in each result directory (from 1 to max K) new subdirectory "admixture" (creating subdirs in Windows in R is complicated)
	system(paste("mkdir ", hauss.geneland.path.mcmc.adm))
	HZ(coordinates=hauss.geneland.coord.utm, geno.dip.codom=hauss.geneland.data, path.mcmc.noadm=hauss.geneland.path.mcmc, nit=hauss.geneland.nit, thinning=hauss.geneland.thinning, path.mcmc.adm=hauss.geneland.path.mcmc.adm)
	}
# Produce maps of respective inferred clusters
for (hauss.geneland.irun in 1:hauss.geneland.nrun) {
	hauss.geneland.path.mcmc <- paste("geneland/", hauss.geneland.irun, "/", sep="")
	# Maps - tesselations
	PlotTessellation(coordinates=hauss.geneland.coord.utm, path.mcmc=hauss.geneland.path.mcmc, printit=TRUE, path=hauss.geneland.path.mcmc)
	for (hauss.geneland.irun.img in 1:hauss.geneland.maxpop) { dev.off() } # We must close graphical device manually
	}
# Estimate frequencies of null alleles
hauss.geneland.fna <- list()
for (hauss.geneland.irun in 1:hauss.geneland.nrun) {
	hauss.geneland.path.mcmc <- paste("geneland/", hauss.geneland.irun, "/", sep="")
	# Estimation
	hauss.geneland.fna[[hauss.geneland.irun]] <- EstimateFreqNA(path.mcmc=hauss.geneland.path.mcmc)
	}
# See output
hauss.geneland.fna
# Calculate average posterior probability
hauss.geneland.lpd <- rep(NA, hauss.geneland.nrun)
for (hauss.geneland.irun in 1:hauss.geneland.nrun) {
	hauss.geneland.path.mcmc <- paste("geneland/", hauss.geneland.irun, "/", sep="")
	hauss.geneland.path.lpd <- paste(hauss.geneland.path.mcmc,"log.posterior.density.txt",sep="")
	hauss.geneland.lpd[hauss.geneland.irun] <- mean(scan(hauss.geneland.path.lpd)[-(1:hauss.geneland.burnin)])
	}
# Sorts runs according to decreasing posterior probability - the first one is the best
order(hauss.geneland.lpd, decreasing=TRUE)
hauss.geneland.lpd

## Maps

# Load libraries
library(sp)
library(rworldmap) # Basic world maps
library(TeachingDemos) # To be able to move text little bit
library(RgoogleMaps) # Google maps
library(mapplots) # Plot pie charts
library(adegenet)
# Plot basic map with state boundaries within selected range
plot(x=getMap(resolution="high"), xlim=c(19, 24), ylim=c(39, 44), asp=1, lwd=1.5)
box() # Add frame around the map
# Plot location points
points(x=hauss.genpop@other$xy[["lon"]], y=hauss.genpop@other$xy[["lat"]], pch=15:19, col="red", cex=4)
# Add text descriptions for points. Text is aside and with background
shadowtext(x=hauss.genpop@other$xy[["lon"]], y=hauss.genpop@other$xy[["lat"]], labels=as.vector(popNames(hauss.genind)), col="black", bg="white", theta=seq(pi/4, 2*pi, length.out=8), r=0.15, pos=c(1, 3, 2, 4, 4), offset=0.75, cex=1.5)
# Insert legend
legend(x="topright", inset=1/50, legend=c("He", "Oh", "Pr", "Ne", "Sk"), col="red", border="black", pch=15:19, pt.cex=2, bty="o", bg="lightgrey", box.lwd=1.5, cex=1.5, title="Populations")

# Google map is produced into a file
# Google recently started to require API key, see https://developers.google.com/maps/documentation/maps-static/overview
?GetMap # See all options
?PlotOnStaticMap # See all options
# Download map
hauss.gmap <- GetMap(center=c(lat=41, lon=21), size=c(640, 640), destfile="gmap.png", zoom=8, maptype="satellite")
# Plot saved map, with extra data
PlotOnStaticMap(MyMap=hauss.gmap, lat=hauss.genpop@other$xy[["lat"]], lon=hauss.genpop@other$xy[["lon"]], FUN=points, pch=19, col="blue", cex=5)
PlotOnStaticMap(MyMap=hauss.gmap, lat=hauss.genpop@other$xy[["lat"]], lon=hauss.genpop@other$xy[["lon"]], add=TRUE, FUN=points, pch=19, col="red", cex=3)
PlotOnStaticMap(MyMap=hauss.gmap, lat=hauss.genpop@other$xy[["lat"]], lon=hauss.genpop@other$xy[["lon"]], add=TRUE, FUN=text, labels=as.vector(popNames(hauss.genind)), pos=4, cex=3, col="white")
# Google maps have their own internal scaling, adding of points by standard functions will not work correctly

# Pie charts on a map
# Prepare matrix with some data (exemplary distribution of haplotypes)
hauss.pie <- cbind(c(20, 30, 15, 40, 10), c(30, 10, 25, 5, 45), c(10, 20, 20, 15, 5))
rownames(hauss.pie) <- popNames(hauss.genpop) # Add row names according to populations
colnames(hauss.pie) <- c("HapA", "HapB", "HapC") # Add column names according to data displayed
class(hauss.pie) # Check it is matrix
hauss.pie # See resulting matrix
# Plot basic map with state boundaries within selected range
plot(x=getMap(resolution="high"), xlim=c(20, 23), ylim=c(41, 42), asp=1, lwd=1.5)
box() # Add frame around the map
# Plot the pie charts
for (L in 1:5) { add.pie(z=hauss.pie[L,], x=as.vector(hauss.genpop@other$xy[["lon"]])[L], y=as.vector(hauss.genpop@other$xy[["lat"]])[L], labels=names(hauss.pie[L,]), radius=0.1, col=topo.colors(3)) }
?add.pie # See more options
# Add population labels
text(x=hauss.genpop@other$xy[["lon"]], y=hauss.genpop@other$xy[["lat"]], labels=as.vector(popNames(hauss.genind)), col="red", cex=2)

# Pie charts on Google map
# Prepare list to store recalculated coordinates
hauss.gmap.coord <- list()
# Calcualtion of coordinates to form required by Google Maps
for (LC in 1:5) { hauss.gmap.coord[[LC]] <- LatLon2XY.centered(MyMap=hauss.gmap, lat=as.vector(hauss.genpop@other$xy[["lat"]])[LC], lon=as.vector(hauss.genpop@other$xy[["lon"]])[LC], zoom=8) }
hauss.gmap.coord # See result
# Plot plain map
PlotOnStaticMap(MyMap=hauss.gmap)
# Plot pie charts
for (LP in 1:5) { add.pie(z=hauss.pie[LP,], x=hauss.gmap.coord[[LP]]$newX, y=hauss.gmap.coord[[LP]]$newY, labels=names(hauss.pie[LP,]), radius=25, col=topo.colors(n=3, alpha=0.7)) }
# Alternative option to plot pie charts
for (LF in 1:5) { plotrix::floating.pie(xpos=hauss.gmap.coord[[LF]]$newX, ypos=hauss.gmap.coord[[LF]]$newY, x=hauss.pie[LF,], radius=30, col=heat.colors(n=3, alpha=0.5) ) }
# Add population text labels
PlotOnStaticMap(MyMap=hauss.gmap, lat=hauss.genpop@other$xy[["lat"]], lon=hauss.genpop@other$xy[["lon"]], add=TRUE, FUN=text, labels=as.vector(popNames(hauss.genind)), cex=2.5, col="white")

library(maps) # Various mapping tools (plotting, ...)
library(mapdata) # More detailed maps, but political boundaries often outdated, see https://CRAN.R-project.org/package=mapdata
library(mapproj) # Converts latitude/longitude into projected coordinates
# Plot a map, check parameters, among others projection and ?mapproject for its details
map(database="world", boundary=TRUE, interior=TRUE, fill=TRUE, col="lightgrey", plot=TRUE, xlim=c(16, 27), ylim=c(37, 46))
# If you'd use projection, use - mapproject() to convert also coordinates! See ?mapproject for details
points(x=hauss.genpop@other$xy[["lon"]], y=hauss.genpop@other$xy[["lat"]], pch=15:19, col="red", cex=3)

# Plotting on SHP files
library(maptools)
library(rgdal)
# Load SHP file
# Data from https://download.geofabrik.de/europe/macedonia.html
# R working directory has to contain also respective DBF and SHX files (same name, only different suffix)
dir() # Verify required files are unpacked in the working directory
# Get from https://soubory.trapa.cz/rcourse/macedonia.zip
# Check correct import by plotting all layers
macedonia_building <- readOGR(dsn="macedonia_buildings.shp")
plot(macedonia_building)
macedonia_landuse <- readOGR(dsn="macedonia_landuse.shp")
plot(macedonia_landuse)
macedonia_natural <- readOGR(dsn="macedonia_natural.shp")
plot(macedonia_natural)
macedonia_railways <- readOGR(dsn="macedonia_railways.shp")
plot(macedonia_railways)
macedonia_roads <- readOGR(dsn="macedonia_roads.shp")
plot(macedonia_roads)
macedonia_waterways <- readOGR(dsn="macedonia_waterways.shp")
plot(macedonia_waterways)
# Plot all layers into single image, add more information
plot(macedonia_building)
plot(macedonia_landuse, add=TRUE, col="darkgreen", fill=TRUE)
plot(macedonia_natural, add=TRUE, col="green", fill=TRUE)
plot(macedonia_railways, add=TRUE, col="brown", lty="dotted")
plot(macedonia_roads, add=TRUE, col="orange")
plot(macedonia_waterways, add=TRUE, col="blue", lwd=2)
# Add state boundaries
plot(x=getMap(resolution="high"), xlim=c(19, 24), ylim=c(39, 44), asp=1, lwd=5, add=TRUE) # Or e.g.
map(database="world", boundary=TRUE, interior=TRUE, fill=FALSE, col="red", add=TRUE, plot=TRUE, xlim=c(16, 27), ylim=c(37, 46), lwd=5)
# Add sampling points
points(x=hauss.genpop@other$xy[["lon"]], y=hauss.genpop@other$xy[["lat"]], pch=15:19, col="red", cex=4)
# Add description of psampling points
shadowtext(x=hauss.genpop@other$xy[["lon"]], y=hauss.genpop@other$xy[["lat"]], labels=as.vector(popNames(hauss.genind)), col="black", bg="white", theta=seq(pi/4, 2*pi, length.out=8), r=0.15, pos=c(1, 3, 2, 4, 4), offset=0.75, cex=1.5)
# Add legend
legend(x="topright", inset=1/50, legend=c("He", "Oh", "Pr", "Ne", "Sk"), col="red", border="black", pch=15:19, pt.cex=2, bty="o", bg="lightgrey", box.lwd=1.5, cex=1.5, title="Populations")

## Tree manipulations

# Library
library(ape)

# Root the tree
plot.phylo(hauss.nj)
print.phylo(hauss.nj)
tiplabels() # Shows tip numbers
# root.phylo accepts either tip number or tip label
# resolve.root=TRUE ensures root will be bifurcating (without this parameter it sometimes doesn't work)
# outgroup can be single value or vector of multiple tips
hauss.nj.rooted <- root.phylo(phy=hauss.nj, resolve.root=TRUE, outgroup=16) # Or
hauss.nj.rooted <- root.phylo(phy=hauss.nj, resolve.root=TRUE, outgroup="H16") # Or
hauss.nj.rooted <- root.phylo(phy=hauss.nj, resolve.root=TRUE, outgroup=c("H42", "H43"))
# root.multiPhylo() is an alias for root.phylo() for class multiPhylo
print.phylo(hauss.nj.rooted)
plot.phylo(hauss.nj.rooted)
# Check if it is rooted - returns TRUE/FALSE
is.rooted(hauss.nj.rooted)
is.rooted(hauss.nj)
# Root the tree interactively
plot.phylo(hauss.nj)
hauss.nj.rooted <- root.phylo(phy=hauss.nj, interactive=TRUE) # Click to selected tip
plot.phylo(hauss.nj.rooted)
# Unroot the tree
?unroot.phylo
# root(), unroot() and is.rooted() works with single or multiple trees (class phylo or multiPhylo)

# Rotate (swap) clade
# plot.phylo plots tree in exact order as it is in the phylo object
plot.phylo(hauss.nj)
nodelabels()
hauss.nj.rotated <- ape::rotate(phy=hauss.nj, node=74)
plot.phylo(hauss.nj.rotated)
nodelabels()

# Ladderize the tree - step-wise rotation of nodes
plot.phylo(hauss.nj)
hauss.nj.ladderized <- ladderize(hauss.nj)
plot.phylo(hauss.nj.ladderized)

# Drop "extinct" tips - those who don't reach end the tree
# tolerance is respective to the used metrics
plot.phylo(hauss.nj)
axisPhylo()
hauss.nj.fossil <- drop.fossil(phy=hauss.nj, tol=0.4)
plot.phylo(hauss.nj.fossil)
# See details
?drop.fossil

# Extract clade, part of tree
# Plot source tree
plot.phylo(hauss.nj)
# See node labels (numbers) - needed for some tasks
nodelabels()
# Non-interactively extract clade
hauss.nj.extracted <- extract.clade(phy=hauss.nj, node=60)
# See new extracted tree
plot.phylo(hauss.nj.extracted)
# Interactive extraction
plot.phylo(hauss.nj)
# Select clade to extract by clicking on it
hauss.nj.extracted <- extract.clade(phy=hauss.nj, interactive=TRUE)
# See new extracted tree
plot.phylo(hauss.nj.extracted)
# See more options
?extract.clade
?keep.tip # Does the opposite - keeps only selected tip(s)

# Bind two trees into one
# Bind donor "y" tree to a given position of the "x" tree
?bind.tree # See options
plot.phylo(hauss.nj.fossil)
nodelabels()
plot.phylo(hauss.nj.extracted)
nodelabels()
hauss.nj.bind <- bind.tree(x=hauss.nj.fossil, y=hauss.nj.extracted, where="root", position=0, interactive=FALSE)
plot.phylo(hauss.nj.bind)
# Bind two trees interactively
# Plot tree receiving the new one
plot.phylo(hauss.nj.fossil)
# Select where to bind new tree to
hauss.nj.bind <- bind.tree(x=hauss.nj.fossil, y=hauss.nj.extracted, interactive=TRUE)
plot.phylo(hauss.nj.bind)

# Work with multiple trees
# Read tree(s) in NEWICK format - single or multiple tree(s)
oxalis.trees <- read.tree(file="https://soubory.trapa.cz/rcourse/oxalis.nwk")
oxalis.trees
lapply(X=oxalis.trees, FUN=print.phylo)
plot.multiPhylo(x=oxalis.trees) # See all trees
summary(oxalis.trees)
length(oxalis.trees)
names(oxalis.trees)
# Export trees in NEWICK format
write.tree(phy=oxalis.trees, file="trees.nwk")
# Export trees in NEXUS format
write.nexus(oxalis.trees, file="trees.nexus")
# Root all trees
oxalis.trees.rooted <- root.multiPhylo(phy=oxalis.trees, outgroup="O._fibrosa_S159", resolve.root=TRUE)
lapply(X=oxalis.trees.rooted, FUN=print.phylo)
plot.multiPhylo(x=oxalis.trees) # See all trees

# Drop a tip
plot.phylo(hauss.nj)
hauss.nj[["tip.label"]]
tiplabels()
hauss.nj.drop <- drop.tip(phy=hauss.nj, tip=47) # Or
hauss.nj.drop <- drop.tip(phy=hauss.nj, tip="H31") # Or
hauss.nj.drop <- drop.tip(phy=hauss.nj, tip=c("H18", "H29", "H31"))
plot.phylo(hauss.nj.drop)

# Drop a tip from multiPhylo
plot.multiPhylo(x=oxalis.trees)
# See tip labels
oxalis.trees[[1]][["tip.label"]]
oxalis.trees.drop <- lapply(X=oxalis.trees, FUN=drop.tip, tip="O._callosa_S15")
class(oxalis.trees.drop) <- "multiPhylo"
plot.multiPhylo(x=oxalis.trees.drop)
lapply(X=oxalis.trees.drop, FUN=print.phylo)

# Check if the tree is ultrametric - is variance of distances of all tips to node 0? It is required for some analysis
is.ultrametric(oxalis.trees)

# Fitting a chronogram to a phylogenetic tree whose branch lengths are in number of substitution per sites - force tree to be ultrametric
?chronos

# Compute branch lengths for trees without branch lengths
?compute.brlen

# Computes the branch lengths of a tree giving its branching times (aka node ages or heights)
?compute.brtime

## Maximum parsimony
library(phangorn)
# Conversion to phyDat for phangorn
gunnera.phydat <- as.phyDat(gunnera.mafft.ng)
# Prepare starting tree
gunnera.tre.ini <- nj(dist.dna(x=gunnera.mafft.ng, model="raw"))
# Parsimony details
?parsimony
# Returns maximum parsimony score
parsimony(tree=gunnera.tre.ini, data=gunnera.phydat)
# Optimisation - returns maximum parsimony tree
gunnera.tre.pars <- optim.parsimony(tree=gunnera.tre.ini, data=gunnera.phydat)
# Draw a tree
plot.phylo(x=gunnera.tre.pars, type="cladogram", edge.width=2)
title("Maximum-parsimony tree of Gunnera spp.")

## Topographical distances among trees

library(gplots)
library(corrplot)
library(phytools)

# Compute matrix of topological distances among phylogenetic trees
?dist.topo # See details of available computing methods
oxalis.trees.d <- dist.topo(x=oxalis.trees, method="score")

# Basic information about the distance matrix
dim(as.matrix(oxalis.trees.d))
head.matrix(as.matrix(oxalis.trees.d))

# Create heatmaps using heatmap.2 function from gplots package
heatmap.2(x=as.matrix(x=oxalis.trees.d), Rowv=FALSE, Colv="Rowv", dendrogram="none", symm=TRUE, scale="none", na.rm=TRUE, revC=FALSE, col=rainbow(15), cellnote=round(x=as.matrix(x=oxalis.trees.d), digits=2), notecex=1, notecol="white", trace="row", linecol="black", labRow=names(oxalis.trees), labCol=names(oxalis.trees), key=TRUE, keysize=2, density.info="density", symkey=FALSE, main="Correlation matrix of topographical distances", xlab="Trees", ylab="Trees")

# Robinsons-Foulds distance
oxalis.trees.d.rf <- multiRF(oxalis.trees)

# Add names of columns and rows
colnames(oxalis.trees.d.rf) <- names(oxalis.trees)
rownames(oxalis.trees.d.rf) <- names(oxalis.trees)

# Create heatmap using corrplot function from corrplot package
corrplot(corr=oxalis.trees.d.rf, method="circle", type="upper", col=rainbow(15), title="Correlation matrix of topographical distances", is.corr=FALSE, diag=FALSE, outline=TRUE, order="alphabet", tl.pos="lt", tl.col="black")
corrplot(corr=oxalis.trees.d.rf, method="number", type="lower", add=TRUE, col=rainbow(15), title="Correlation matrix of topographical distances", is.corr=FALSE, diag=FALSE, outline=FALSE, order="alphabet", tl.pos="ld", cl.pos="n")

# PCoA from distance matrices of topographical differences among trees
# Test if the distance matrix is Euclidean or not
library(ade4)
is.euclid(distmat=oxalis.trees.d, plot=TRUE)
# Calculate the PCoA
oxalis.trees.pcoa <- dudi.pco(d=oxalis.trees.d, scannf=TRUE, full=TRUE)
# Plot PCoA
s.label(dfxy=oxalis.trees.pcoa$li)
# Add kernel densities
s.kde2d(dfxy=oxalis.trees.pcoa$li, cpoint=0, add.plot=TRUE)
# Add histogram of eigenvalues
add.scatter.eig(oxalis.trees.pcoa[["eig"]], 3,1,2, posi="topleft")
# Add title to the plot
title("PCoA of matrix of pairwise trees distances")
# Alternative function to plot PCA plot
scatter(x=oxalis.trees.pcoa, posieig="topleft")
# See current state
oxalis.trees
# Remove outlying trees
oxalis.trees[c("T471", "T639", "T654")] <- NULL
# See after removal
oxalis.trees

## kdetrees
# Load library
library(kdetrees)
?kdetrees # See options
# Run main function - play with parameter k
oxalis.kde <- kdetrees(trees=oxalis.trees, k=0.4, distance="dissimilarity", topo.only=FALSE, greedy=TRUE)
# See results
oxalis.kde
plot(oxalis.kde)
hist(oxalis.kde)
# See removed trees
plot.multiPhylo(oxalis.kde[["outliers"]])
# Save removed trees
write.tree(phy=oxalis.kde[["outliers"]], file="oxalis_trees_outliers.nwk")
# Save kdetrees report
write.table(x=as.data.frame(x=oxalis.kde), file="oxalis_trees_scores.tsv", quote=FALSE, sep="\t")
# Extract passing trees
oxalis.trees.good <- oxalis.trees[names(oxalis.trees) %in% names(oxalis.kde[["outliers"]]) == FALSE]
oxalis.trees.good
# Save passing trees
write.tree(phy=oxalis.trees.good, file="trees_good.nwk")

## Seeing trees in the forest

# Consenus tree (50% rule)
oxalis.tree.con <- ape::consensus(oxalis.trees.rooted, p=0.5, check.labels=TRUE)
print.phylo(oxalis.tree.con)
# Plot the tree
plot.phylo(oxalis.tree.con, edge.width=2, label.offset=0.3)
axisPhylo(side=1)

# Species tree
# All trees must be ultrametric - chronos scale them
oxalis.trees.ultra <- lapply(X=oxalis.trees.rooted, FUN=chronos, model="correlated")
class(oxalis.trees.ultra) <- "multiPhylo"
# Mean distances
?speciesTree # See underlying method
oxalis.tree.sp.mean <- speciesTree(x=oxalis.trees.ultra, FUN=mean)
print.phylo(oxalis.tree.sp.mean)
# Plot the tree
plot.phylo(oxalis.tree.sp.mean, edge.width=2, label.offset=0.01)
edgelabels(text=round(oxalis.tree.sp.mean[["edge.length"]], digits=2), frame="none", col="red", bg="none")
axisPhylo(side=1)

# Parsimony super tree
library(phytools)
?mrp.supertree # See details
oxalis.tree.sp <- mrp.supertree(tree=oxalis.trees.rooted, method="optim.parsimony", rooted=TRUE)
print.phylo(oxalis.tree.sp)
plot.phylo(oxalis.tree.sp, edge.width=2, label.offset=0.01)
axisPhylo(side=1)
?phangorn::superTree # Similar function
?phangorn::coalSpeciesTree # Has coalescence model to handle multiple individuals per species

# TODO STAR and STEAC species trees
# https://github.com/lliu1871/phybase

# Networks
library(phangorn)
oxalis.tree.net <- consensusNet(oxalis.trees.rooted, prob=0.25)
plot(x=oxalis.tree.net, planar=FALSE, type="2D", use.edge.length=TRUE, show.tip.label=TRUE, show.edge.label=TRUE, show.node.label=TRUE, show.nodes=TRUE, edge.color="black", tip.color="blue") # 2D
plot(x=oxalis.tree.net, planar=FALSE, type="3D", use.edge.length=TRUE, show.tip.label=TRUE, show.edge.label=TRUE, show.node.label=TRUE, show.nodes=TRUE, edge.color="black", tip.color="blue") # 3D

# Density tree
# Prepare list of trees to show
hauss.nj.trees <- list(hauss.nj, hauss.nj.bruvo, hauss.nj.rooted)
hauss.nj.trees <- lapply(X=hauss.nj.trees, FUN=compute.brlen)
hauss.nj.trees <- lapply(X=hauss.nj.trees, FUN=chronos)
class(hauss.nj.trees) <- "multiPhylo"
# The trees should be (otherwise plotting works, but may be more ugly)...
is.rooted.multiPhylo(hauss.nj.trees) # rooted,
is.ultrametric.multiPhylo(hauss.nj.trees) # ultrametric and
is.binary.multiPhylo(hauss.nj.trees) # binary bifurcating.
# Plotting has various options, play with it
phangorn::densiTree(x=hauss.nj.trees, direction="downwards", scaleX=TRUE, col=rainbow(3), width=5, cex=1.5)
densiTree(x=hauss.nj.trees, direction="upwards", scaleX=TRUE, width=5)
densiTree(x=hauss.nj.trees, scaleX=TRUE, width=5, cex=1.5)
?phangorn::densiTree
?phytools::densityTree
# Another version for comparing several trees
phytools::densityTree(trees=oxalis.trees.ultra, fix.depth=TRUE, use.gradient=TRUE, alpha=0.5, lwd=4) # Probably to much noise... :-?
phytools::densityTree(trees=oxalis.trees.ultra[1:3], fix.depth=TRUE, use.gradient=TRUE, alpha=0.5, lwd=4) # Nice selection
phytools::densityTree(trees=oxalis.trees.ultra[c(2,4,6,7)], fix.depth=TRUE, use.gradient=TRUE, alpha=0.5, lwd=4) # Nice selection

# Plot all trees on same scale
kronoviz(x=oxalis.trees.rooted, layout=length(oxalis.trees.rooted), horiz=TRUE)
# Close graphical device to cancel division of plotting device
dev.off()

## Compare two trees
# Compare topology of the species trees - basically outputs TRUE/FALSE
all.equal.phylo(oxalis.tree.sp, oxalis.tree.sp.mean, use.edge.length=FALSE)
?all.equal.phylo # Use to see comparison possibilities
# Plot two trees with connecting lines
# We need 2 column matrix with tip labels
tips.labels <- matrix(data=c(sort(oxalis.tree.sp[["tip.label"]]), sort(oxalis.tree.sp.mean[["tip.label"]])), nrow=length(oxalis.tree.sp[["tip.label"]]), ncol=2)
tips.labels
# Draw a tree - play with graphical parameters and use rotate=TRUE
# to be able to adjust fit manually
cophyloplot(x=ladderize(oxalis.tree.sp), y=ladderize(oxalis.tree.sp.mean), assoc=tips.labels, use.edge.length=FALSE, space=60, length.line=1, gap=2, type="phylogram", rotate=TRUE, col="red", lwd=1.5, lty=2)
title("Comparing the trees\nParsimony super tree\tSpecies tree")
legend("topleft", legend="Red lines\nconnect tips", text.col="red", cex=0.75, bty="n", x.intersp=-2, y.intersp=-2)
# Alternative implementation
oxalis.cophylo <- cophylo(tr1=oxalis.tree.sp, tr2=oxalis.tree.sp.mean, assoc=(cbind(sort(oxalis.tree.sp$tip.label), sort(oxalis.tree.sp$tip.label))), rotate=TRUE)
plot.cophylo(x=oxalis.cophylo, lwd=2, link.type="curved")
title("Comparisson of species tree (left) and parsimony supertree (right)")

## More about plotting the trees
?plot.phylo # check it for various possibilities what to influence
?par
?points
par(mfrow=c(1, 2)) # Plot two plots in one row
plot.phylo(x=hauss.nj, type="cladogram", use.edge.length=FALSE, direction="rightwards")
plot.phylo(x=hauss.nj, type="cladogram", use.edge.length=FALSE, direction="leftwards")
dev.off() # Close graphical device to cancel par() settings

# Nice tiplabels and higlighted tiplabel
# Load tree in text format
trape <- read.tree(text="((Homo, Pan), Gorilla);")
# Plot the tree
plot.phylo(x=trape, show.tip.label=FALSE)
# Add colored tip labels
tiplabels(trape[["tip.label"]], bg=c("white", "black", "white"), col=c("black", "white", "black"), cex=2)
# Add colored node labels
nodelabels(text=c("6.4 Ma", "5.4 Ma"), frame="circle", bg="yellow")
# Add scale bar
add.scale.bar()
# Note vectors for tip/node labels

## PIC

# Load library
library(ape)

# Load training data
data(shorebird, package="caper")
# See it
?caper::shorebird
# See the data part
head(shorebird.data)
# See the phylogeny
shorebird.tree
plot.phylo(shorebird.tree)
# The tree must be fully bifurcating for most of methods
# multi2di randomly splits polytomies into bifurcating topology
shorebird.tree <- multi2di(phy=shorebird.tree)
# See result
plot.phylo(shorebird.tree)

# PIC - phylogenetically independent contrasts
shorebird.pic.fmass <- pic(x=shorebird.data[["F.Mass"]], phy=shorebird.tree, scaled=TRUE, var.contrasts=FALSE, rescaled.tree=FALSE)
shorebird.pic.eggmass <- pic(x=shorebird.data[["Egg.Mass"]], phy=shorebird.tree, scaled=TRUE, var.contrasts=FALSE, rescaled.tree=FALSE)

# Plot a tree with PIC values
plot.phylo(x=shorebird.tree, edge.width=2, cex=0.8)
nodelabels(round(shorebird.pic.fmass, digits=3), adj=c(0, -0.5), frame="none", col="red")
nodelabels(round(shorebird.pic.eggmass, digits=3), adj=c(0, 1), frame="none", col="red")
add.scale.bar()

# Plot PIC
plot(x=shorebird.pic.fmass, y=shorebird.pic.eggmass, pch=16, cex=1.5)
abline(a=0, b=1, lty=2) # x=y line

# Correlation of PIC of body mass and longevity
cor(x=shorebird.pic.fmass, y=shorebird.pic.eggmass, method="pearson")
# Correlation test
cor.test(x=shorebird.pic.fmass, y=shorebird.pic.eggmass, alternative="greater", method="pearson")
# Linear model of both PICs
lm(formula=shorebird.pic.fmass~shorebird.pic.eggmass)
# Because PICs have expected mean zero - such linear regressions should be done through the origin (i.e. the intercept is set to zero)
lm(formula=shorebird.pic.fmass~shorebird.pic.eggmass-1)

# Permutation procedure to test PIC
lmorigin(formula=shorebird.pic.fmass~shorebird.pic.eggmass, nperm=1000)

# Intraspecific variation
# PIC - orthonormal contrasts
# List of numerical vectors is required as an input, otherwise usage is same as with pic()
?pic.ortho

## Phylogenetic autocorrelation

# Autocorrelation coefficient to quantify whether the distribution of a trait among a set of species is affected or not by their phylogenetic relationships
# In the absence of phylogenetic autocorrelation, the mean expected value of I and its variance are known - it is thus possible to test the null hypothesis of the absence of dependence among observations
# Let's choose weights as wij = 1/dij, where the d’s is the distances measured on the tree - cophenetic() calculates cophenetic distances
shorebird.weights <- 1/cophenetic(shorebird.tree) # can be just cophenetic(shorebird.tree) or some other transformation
# See it
class(shorebird.weights)
head(shorebird.weights)
# Set diagonal to 0
diag(shorebird.weights) <- 0
# Calculate Moran's I
Moran.I(x=shorebird.data[["M.Mass"]], weight=shorebird.weights, alternative="greater") # Significant, but super small, positive phylogenetic correlation among body mass
Moran.I(x=shorebird.data[["F.Mass"]], weight=shorebird.weights, alternative="greater") # Significant, but super small, positive phylogenetic correlation among body mass

# Test of Moran's with randomisation procedure
library(ade4)
?gearymoran
gearymoran(bilis=shorebird.weights, X=shorebird.data[,2:5], nrepet=1000) # For all characters significant, but very small

# Test of Abouheif designed to detect phylogenetic autocorrelation in a quantitative trait - in fact Moran's I test using a particular phylogenetic proximity between tips
library(adephylo)
abouheif.moran(x=shorebird.data[,2:5], W=shorebird.weights, method="oriAbouheif", nrepet=1000, alter="greater")

# correlogram can be used to visualize the results of phylogenetic autocorrelative analysis
# Loads training data set
data(carnivora)
# Look at the data
head(carnivora)
# Calculate the correlogram
?correlogram.formula
carnivora.correlogram <- correlogram.formula(formula=SW~Order/SuperFamily/Family/Genus, data=carnivora)
# See results
carnivora.correlogram
# Calculate the correlogram - test for both body masses
carnivora.correlogram2 <- correlogram.formula(formula=SW+FW~Order/SuperFamily/Family/Genus, data=carnivora)
# See results
carnivora.correlogram2
# Plot it
plot(x=carnivora.correlogram, legend=TRUE, test.level=0.05, col=c("white", "black"))
# Plot it - test for both body masses
# Two graphs
plot(x=carnivora.correlogram2, lattice=TRUE, legend=TRUE, test.level=0.05)
# Only one graph
plot(x=carnivora.correlogram2, lattice=FALSE, legend=TRUE, test.level=0.05)

## Phylogenetic PCA
library(adephylo) # Library needed to create phylo4d object required by ppca
# Calculate pPCA
shorebird.ppca <- ppca(x=phylo4d(x=shorebird.tree, shorebird.data[,2:5]), method="patristic", center=TRUE, scale=TRUE, scannf=TRUE, nfposi=1, nfnega=0)
# Print results
print(shorebird.ppca)
# See summary information
summary(shorebird.ppca)
# See PCA scores for variables on phylogenetic tree
scatter(shorebird.ppca)
# See decomposition of pPCA eigenvalues
screeplot(shorebird.ppca)
# Plot pPCA results - global vs. local structure, decomposition of pPCA eigenvalues, PCA plot of variables and PCA scores for variables on phylogenetic tree
plot(shorebird.ppca)

## Orthonormal Decomposition - phylogenetic eigenvector regression

# Decomposition of topographical distances (right plot)
library(phylobase)
table.phylo4d(x=phylo4d(x=shorebird.tree, tip.data=treePart(x=shorebird.tree, result="orthobasis")), treetype="cladogram")

# Generate some random variable
library(geiger)
shorebird.eco <- sim.char(phy=shorebird.tree, par=0.1, nsim=1, model="BM")[,,1]
?sim.char # See it for another possibilities to simulate data
# Names for the values
names(shorebird.eco) <- shorebird.tree[["tip.label"]]
shorebird.eco # See it

# significant result - significant phylogenetic inertia (phylogenetic effect) - the tendency for traits to resist evolutionary change despite environmental perturbations
anova(lm(shorebird.eco ~ as.matrix(orthobasis.phylo(x=shorebird.tree, method="patristic")[,1:2])))
anova(lm(shorebird.data[["M.Mass"]] ~ as.matrix(orthobasis.phylo(x=shorebird.tree, method="patristic")[,1:2])))
orthogram(x=shorebird.eco, tre=shorebird.tree, nrepet=1000, alter="two-sided")
?orthogram # See another calculation possibilities
orthogram(x=shorebird.data[["M.Mass"]], tre=shorebird.tree, nrepet=1000, alter="two-sided")

## Phylogenetic Generalized Least Squares

# Functions corBlomberg, corBrownian, corMartins and corPagel from ape package create correlation matrix of evolution of continuous character according to the given tree - same parameters (except 'value')
# Fitting an Ornstein-Uhlenbeck Motion model in PGLS
library(nlme)
summary(gls(model=F.Mass ~ Egg.Mass, data=shorebird.data, correlation=corBrownian(value=1, phy=shorebird.tree)))

# Implementation in caper package
library(caper) # Load needed library
shorebird.pgls <- pgls(formula=shorebird.data[["F.Mass"]] ~ shorebird.data[["Egg.Mass"]], data=comparative.data(phy=shorebird.tree, data=as.data.frame(cbind(shorebird.data[["F.Mass"]], shorebird.data[["Egg.Mass"]], shorebird.data[["Species"]])), names.col=V3, vcv=TRUE))
summary(shorebird.pgls) # See the result
# See the plot of observer and fitted values
plot(shorebird.pgls)
abline(a=0, b=1, col="red")
anova(shorebird.pgls) # ANOVA view of the model
AIC(shorebird.pgls) # Akaike's information criterion (smaller = better)

## Generalized Estimating Equations

# Calculate the model
compar.gee(formula=shorebird.data[["F.Mass"]] ~ shorebird.data[["Egg.Mass"]], phy=shorebird.tree)
# or with correlation matrix
compar.gee(formula=shorebird.data[["F.Mass"]] ~ shorebird.data[["Egg.Mass"]], corStruct=corMartins(value=1, phy=shorebird.tree, fixed=TRUE))
# for corStruct there are similar functions corBlomberg, corMartins, corPagel, corBrownian - see manuals for differences

# # TODO multiple phylogenetic regressions and residuals
# # x can be matrix of data
# phyl.resid(tree=shorebird.tree, x=primates.body, Y=primates.longevity, method="BM")
# phyl.resid(tree=shorebird.tree, x=primates.body, Y=primates.longevity, method="lambda")

# # TODO Ornstein-Uhlenbeck Model
# # Simulation of evolution on phylogenetic tree
# compar.ou(x=shorebird.eco, phy=shorebird.tree, alpha=100)

# ## TODO Variance Partitioning
# # vare: the estimated residual variance–covariance matrix;
# # vara: the estimated additive effect variance–covariance matrix;
# # u: the estimates of the phylogeny wide means;
# # A: the additive value estimates;
# # E: the residual value estimates;
# # lik: the log-likelihood.
# # x can be vector, matric or df
# apiaceae.lynch <- compar.lynch(x=shorebird.eco, G=vcv.phylo(phy=shorebird.tree, model="Brownian", corr=TRUE))
# mantel.test(m1=apiaceae.lynch$vara, m2=apiaceae.lynch$vare, nperm=1000, graph=TRUE, alternative="two.sided")

## Phylogenetic Signal

library(picante)
# If Blomberg's values of 1 correspond to a Brownian motion process, which implies some degree of phylogenetic signal or conservatism. K values closer to zero correspond to a random or convergent pattern of evolution, while K values greater than 1 indicate strong phylogenetic signal and conservatism of traits.
# It requires named vector of trait values
shorebird.mmass <- shorebird.data[["M.Mass"]]
names(shorebird.mmass) <- rownames(shorebird.data)
# Bloomberg's K statistics
Kcalc(x=shorebird.mmass, phy=shorebird.tree, checkdata=TRUE)
# Test with permutations
phylosignal(x=shorebird.mmass, phy=shorebird.tree, reps=1000, checkdata=TRUE)
# Analyze multiple traits in once with multiPhylosignal, it requires data frame of numerical values
multiPhylosignal(x=shorebird.data[,2:5], phy=shorebird.tree, reps=1000)

# When there are vectors with standard errors of measurements
library(phytools)
?phylosig # See for details
# Test for phylogenetic signal (here without SE)
phylosig(tree=shorebird.tree, x=shorebird.mmass, method="K", test=TRUE, nsim=1000)
plot(phylosig(tree=shorebird.tree, x=shorebird.mmass, method="K", test=TRUE, nsim=1000))
phylosig(tree=shorebird.tree, x=shorebird.mmass, method="lambda", test=TRUE)
# phylosig() can be used as an alternative to phylosignal() - the functions are similar in basic usage

# Examples of usage of GLS for testing of phylogenetic signal
summary(gls(model=shorebird.mmass ~ 1, correlation=corBrownian(value=1, phy=shorebird.tree)))
summary(pgls(formula=shorebird.mmass ~ 1, data=comparative.data(phy=shorebird.tree, data=as.data.frame(cbind(shorebird.data[["M.Mass"]], shorebird.data[["Species"]])), names.col=V2, vcv=TRUE)))

## Ancestral state reconstruction

# Loading data
# Ackerly & Donoghue (1998) https://doi.org/10.1086/286208
data(maples, package="adephylo")
?adephylo::maples

# Process the phylogenetic tree
# maples data provide tree as plain text in NEWICK, must be imported into the phylo object
maples.tree <- read.tree(text=maples[["tre"]])
maples.tree
plot.phylo(maples.tree)
# For plenty of analysis it must be fully resolved (bifurcating), rooting and ultrametricity often helps
is.binary.phylo(maples.tree)
is.rooted.phylo(maples.tree)
is.ultrametric.phylo(maples.tree)

# See the character matrix
head(maples[["tab"]])
maples.data <- maples[["tab"]][,1:30]
head(maples.data)
summary(maples.data)

# Maples mature height (m)
maples.height <- maples[["tab"]][["MatHt"]]
names(maples.height) <- rownames(maples[["tab"]])
maples.height

# Maples seed size (mg)
maples.sdsz  <- maples[["tab"]][["SdSz"]]
names(maples.sdsz) <- rownames(maples[["tab"]])
maples.sdsz

# Maples leaf + petiole length (mm)
maples.lfpt <- maples[["tab"]][["LfPt"]]
names(maples.lfpt) <- rownames(maples[["tab"]])
maples.lfpt

# By default performs estimation for continuous characters assuming a Brownian motion model fit by maximum likelihood
library(ape)
?ace # See for possible settings
maples.height.ace <- ace(x=maples.height, phy=maples.tree, type="continuous", method="REML", corStruct=corBrownian(value=1, phy=maples.tree))
# See result - reconstructions are in $ace slot - to be plotted on nodes - 1st column are node numbers
maples.height.ace
# Plot it
plot.phylo(x=maples.tree, lwd=2, cex=0.75)
tiplabels(maples.height, adj=c(1, 0), frame="none", col="blue", cex=0.75)
nodelabels(round(maples.height.ace[["ace"]], digits=2), frame="none", col="red", cex=0.75)
# ACE returns long numbers - truncate them by e.g. 
round(x=..., digits=3) # "x" is vector with ACE values

# Other implementations are available in packages geiger, phangorn, ape, phytools, ...
# Parsimony based method
?ape::MPR
# For continuous characters using Maximum Likelihood
?geiger::fitContinuous
# For continuous characters using Markov Chain Monte Carlo
?geiger::fitContinuousMCMC
# For discrete characters, various models available
?geiger::fitDiscrete
# Marginal reconstruction of the ancestral character states
?phangorn::ancestral.pml

library(phytools)

# ML estimation of a continuous trait
?fastAnc
maples.height.fa <- fastAnc(tree=maples.tree, x=maples.height, vars=FALSE, CI=TRUE)
maples.height.fa
plot.phylo(x=maples.tree, edge.width=2, cex=0.75)
tiplabels(maples.height, adj=c(1, 0), frame="none", col="blue", cex=0.75)
nodelabels(round(maples.height.fa[["ace"]], digits=2), frame="none", col="red", cex=0.75)

# ACE for Brownian evolution with directional trend
?anc.trend
maples.height.ac <- anc.trend(tree=maples.tree, x=maples.height, maxit=1000000)
maples.height.ac

# ACE for Brownian evolution using likelihood
?anc.ML
maples.height.ml <- anc.ML(tree=maples.tree, x=maples.height, maxit=1000000, model="BM")
maples.height.ml

# Bayesian ancestral character estimation
?anc.Bayes
maples.height.bayes <- anc.Bayes(tree=maples.tree, x=maples.height, ngen=1000000) # Use more MCMC generations
maples.height.bayes
# Get end of ancestral states from Bayesian posterior distribution (it should converge to certain values)
tail(maples.height.bayes[["mcmc"]])
maples.height.bayes[["mcmc"]][10001,3:18]
# Get means of ancestral states from Bayesian posterior distribution
colMeans(maples.height.bayes[["mcmc"]][2001:nrow(maples.height.bayes[["mcmc"]]),as.character(1:maples.tree$Nnode+length(maples.tree$tip.label))])
# Plot the ancestral states from posterior distribution (it should converge to certain values)
plot(maples.height.bayes)
# Plot the tree and reconstructed ancestral states
plot.phylo(x=maples.tree, edge.width=2, cex=2)
nodelabels(round(x=maples.height.bayes[["mcmc"]][10001,3:18], digits=2), cex=1.5)

# Continuous map
?contMap
contMap(tree=maples.tree, x=maples.height)
# Change colors with setMap()
maples.contmap <- setMap(x=contMap(tree=maples.tree, x=maples.height), colors=c("white", "black"))
plot(maples.contmap)
# See ?par for more settings

## Phenogram

library(adephylo)
table.phylo4d(x=phylo4d(x=maples.tree, tip.data=maples.data), treetype="cladogram", symbol="circles", scale=FALSE, ratio.tree=0.5)
table.phylo4d(x=phylo4d(x=maples.tree, tip.data=maples.data), treetype="cladogram", symbol="circles", scale=TRUE, ratio.tree=0.5)
phenogram(tree=maples.tree, x=maples.height, ftype="i", colors="red", main="Maples adult heights")
fancyTree(tree=maples.tree, type="phenogram95", x=maples.height, ftype="i", main="95-percentile of Maples adult heights")
# 2 characters on 2 axis
phylomorphospace(tree=shorebird.tree, X=shorebird.data[,2:3], label="horizontal", lwd=2)
# 3D - 3 characters it a rotating cube
phylomorphospace3d(tree=shorebird.tree, X=shorebird.data[,2:4], label=TRUE)
# 3 characters on 2 axis and ancestral state reconstruction for all of them
fancyTree(tree=shorebird.tree, type="scattergram", X=shorebird.data[,2:4], res=500, ftype="i")
# See manuals for more settings
?fancyTree
?phenogram
?phylomorphospace
?phylomorphospace3d
?contMap
?setMap
?par
?plot

# Plotting traits on trees

# See options for plotting functions
?plotTree.wBars # There are more variants available
?dotTree
?plotSimmap
# Plot the trees
# Tip labels with bars with length proportional to character values
plotTree.wBars(tree=shorebird.tree, x=shorebird.mmass, tip.labels=TRUE)
# Tip labels with dots with size proportional to character values
dotTree(tree=shorebird.tree, x=shorebird.mmass, tip.labels=TRUE, type="cladogram")

# ## TODO Disparity-through-time
# disparity(phy=shorebird.tree, data=primates.body)
# dtt(phy=shorebird.tree, data=primates.body, nsim=1000)
# # lineage-through-time plot
# ltt(tree=shorebird.tree, plot=TRUE, drop.extinct=FALSE, log.lineages=TRUE, gamma=TRUE)

## Graphics
# Output figure will be saved to the disk as OutputFile.png
png(filename="OutputFile.png", width=720, height=720, bg="white")
# Here can go any number of functions making plots...
plot(...) # Whatever...
# When using plotting commands, nothing is shown on the screen
# The final plot(s) will be saved by:
dev.off() # Closes graphical device - needed after use of plotting functions png(), svg(), pdf(), ... followed by any function like plot() to write the file(s) to the disk
# filename="OutFiles_%03d.png" # Returns list of files named OutFiles_001.png, OutFiles_002.png, ... Useful for functions returning more graphs.
?png # These functions have various possibilities to set size, whatever.
?svg # Exact possibilities of all 3 functions vary from system to system
?pdf # according to graphical libraries available in the computer.

## Install packages from GitHub
# Needed library
require(devtools)
dev_mode(on=TRUE)
# Install selected package from GitHub (user/project)
install_github("thibautjombart/adegenet")
# when finished go back to normal version
dev_mode(on=FALSE)

## Functions
MyFunction <- function (x, y) {
	# Any commands can be here...
	x + y
	}

# Use as usually:
MyFunction(5, 8)
MyFunction(1, 4)
MyFunction(x=4, y=7)
MF <- MyFunction(9, 15)
MF # See it works

## Loops

# For loop
# Simplest loop - print value of "i" in each step
# "i" is commonly used for various indexing
for (i in 1:5) { print(i) }

# In every step modify value of variable "X" (add 1 to previous value)
X <- 0 # Set initial value
for (i in 1:10) {
	# Any commands can be here...
	print("Loop turn") # Some message for user
	print(i) # Print number of turn - note how it is increasing
	X <- X+i # Rise value of "X" by current value of "i" (see previous line)
	print(paste("Variable value:", X)) # Print current value of "X"
	}
for (i in 10:5) { print(i) } # Can be descending...

# Work on each item of a list object
# Print length of each sequence in gunnera.dna
for (L in 1:length(gunnera.dna)) {
	print(length(gunnera.dna[[L]]))
	}

# While loop - it is done while the condition is valid
# While value of "Q" is < 5 (starting from 0), print it and add 1
Q <- 0
while (Q < 5) { print(Q <- Q+1) }

## If-else branching
XX <- seq(from=-3, to=6.5, by=0.1)
XX
YY <- c()
for (II in 1:length(XX)) {
	if(XX[II] <= 2) { # Executed for XX <= 2
		YY[II] <- XX[II]^2
		} else { # if(XX[II] > 2) # Executed for XX > 2
			YY[II] <- 6-XX[II]
			}
	}
YY
plot(XX, YY) # See the result

# Or (different possibility to get very same result)
# Note "XX" is reused from the previous variant
CC <- function(AA) {
	if(AA <= 2) { # Executed for XX <= 2
		BB <- AA^2
	} else { # Executed for XX > 2
		BB <- 6-AA
		}
	return(BB) # The output value
	}
CC # Previously, "YY" contained values to plot made by the for loop, here "CC" contains function to by used by sapply() when plotting
plot(sapply(XX, CC)) # See the result

################################################################################
# TODO more topics

## RAxML
gunnera.raxml <- raxml(DNAbin=usflu.dna, N="autoFC", exec="~/bin/")
gunnera.raxml <- phyloch::raxml(x=gunnera.nogaps, runs=10, file="gunnera.raxml", path="~/bin/")

## Extra
axisGeo() # package phyloch - adds scale in geological time, has many options

## Polysat - microsattelites and mixed ploidies
library(combinat)
library(polysat)
polysattest <- read.GeneMapper("https://soubory.trapa.cz/rcourse/GeneMapperExample.txt")
# Basic information
summary(polysattest)
Loci(polysattest) # Information about loci
class(polysattest)
# Add information missing in input data file
Description(polysattest) <- "Dataset for the tutorial" # Description of data set
PopNames(polysattest) <- c("PopA", "PopB", "PopC") # Population names
PopInfo(polysattest) <- rep(x=1:3, each=100) # Which individuals belong to which population
PopInfo(polysattest) # Which individuals belong to which population
Usatnts(polysattest) <- c(2, 3, 2) # Ploidies
Usatnts(polysattest) # Get information about ploidy levels for populations
# Show samples
Samples(polysattest)
# Subsetting of populations:
polysattest2 <- Samples(object=polysattest, populations="PopA")
summary(polysattest2)
# View genotypes (for locus loc1) of samples A1-A20
viewGenotypes(polysattest, samples=paste("A", 1:20, sep=""), loci="loc1")
# Show missing data
find.missing.gen(polysattest)
summary(polysattest)
testmat <- meandistance.matrix(polysattest)
pca <- cmdscale(testmat)
mycol <- c("red", "green", "blue")
plot(pca[,1], pca[,2], col=mycol[PopInfo(polysattest)], + main = "PCA with Bruvo distance")
testmat2 <- meandistance.matrix(polysattest, distmetric=Lynch.distance)
pca2 <- cmdscale(testmat2)
plot(pca2[,1], pca2[,2], col=rep(c("red", "green", "blue"), each=100), + main = "PCA with Lynch distance")
polysattest2 <- deleteSamples(polysattest, c("B59", "C30"))
polysattest2 <- deleteLoci(polysattest2, "loc2")

