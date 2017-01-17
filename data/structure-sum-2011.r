### This is a script to summarize the results from several structure runs.
### Version 210109, works for outputs of Structure 2.2 using the admixture, no
### admixture, correlated and uncorrelated allele frequency models, as well as
### for files from the Bioportal of the University of Oslo with an additional
### RANDOMIZE line.
### It works for datasets containing recessive alleles (e.g. AFLP) and
### for datasets including information about populations of origin (but
### not for results from the USEPOPINFO option)
### It does not work for the linkage model.

### Function to make a table of LnP(D) for a series of structure output files

Structure.table <- function(infile, pop=0, locprior=F) {
runtab <- read.table(infile, sep="\t", header=FALSE, row.names=NULL)
runtab <- runtab[sort.list(runtab[ ,1]), ]
runnumbers <- table(runtab[ ,1])
Knb <- length (runnumbers)

lnprob <- matrix(1, 2, 5)
z <- 1

if(scan(file=as.character(runtab[1, 2]), skip= 13, n=1, what="character")== "Run") hd <- 13
if(scan(file=as.character(runtab[1, 2]), skip= 14, n=1, what="character")== "Run") hd <- 14
if(scan(file=as.character(runtab[1, 2]), skip= 15, n=1, what="character")== "Run") hd <- 15

ran <- 0
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(file=as.character(runtab[1, 2]), skip= (hd+7), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(file=as.character(runtab[1, 2]), skip= (hd+8), n=1, what="character")== "RANDOMIZE") ran <- 1

startl <- (hd+23+ran)
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "NO") startl <- (hd+24+ran) 
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "RECESSIVE") {
if(scan(file=as.character(runtab[1, 2]), skip= (hd+7), n=1, what="character")== "NO") startl <- (hd+25+ran) else startl <- (hd+24+ran)}

if (pop == 0) {
for (i in 1:Knb) {
	
	runfiles <- as.vector(runtab[z:(z+runnumbers[i]-1), 2])
	K <- scan(runfiles[1], skip=(hd+3), n=1)
	lnprobab <- matrix(NA, runnumbers[i], 5)

	for (r in 1: runnumbers[i]) {
		lnprobab[r, 1] <- K
		lnprobab[r, 2] <- r
		lnprobab[r, 3] <- as.numeric(scan(runfiles[r], skip=(startl + 3 + 2*K), n=7, what="character")[7])
		minsum <- min(as.numeric(scan(runfiles[r], skip=(startl - 9), n=K)))
		if (minsum == 0) lnprobab[r, 4] <- "yes" else lnprobab[r, 4]<- "no"
		lnprobab[r, 5] <- runfiles[r]
		}
		        
	lnprob <- rbind(lnprob, lnprobab)
	rm(lnprobab)
                               
	z <- z + runnumbers[i]
	rm(runfiles)
	}	 }

if (locprior == T) startl <- startl + 1
 
if (pop > 0) {         
for (i in 1:Knb) {
	
	runfiles <- as.vector(runtab[z:(z+runnumbers[i]-1), 2])
	K <- scan(runfiles[1], skip=(hd+3), n=1)
	lnprobab <- matrix(NA, runnumbers[i], 5)

	for (r in 1: runnumbers[i]) {
		lnprobab[r, 1] <- K
		lnprobab[r, 2] <- r
		lnprobab[r, 3] <- as.numeric(scan(runfiles[r], skip=(startl + 2 + 2*K + pop), n=7, what="character")[7])
		members <- as.matrix(read.table(file=runfiles[r], nrows=pop, skip=(startl - 8))[ ,2:(K+1)])
		minsum <- min(apply(members, 2, sum))
		if (minsum == 0) lnprobab[r, 4] <- "yes" else lnprobab[r, 4]<- "no"
		lnprobab[r, 5] <- runfiles[r]
		}

	lnprob <- rbind(lnprob, lnprobab)
	rm(lnprobab)

	z <- z + runnumbers[i]
	rm(runfiles)
	}	 }	

xy <- dim (lnprob)
lnprob <- lnprob[3:xy[1], ]
plot(lnprob[ ,1], lnprob[ ,3], xlab="K", ylab="Ln P(D)", axes=FALSE, cex=1.8, frame=TRUE)
axis(1, at=c(1:K))
axis(2)
lnprob <- rbind(c("K", "run", "Ln P(D)", "empty groups", "file name"), lnprob)
rm(xy)
write.table(lnprob, file="lnPData.txt", quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)		
}


### Functions which calculate the coefficient of similarity 
### among a series of structure runs with the same K


# function which orders the result matrices of 2 structure runs 
# the same groups will be in the same columns

ordermat <- function(mat1, mat2, indnr, Ka) {

sr <- rank(apply(mat2, 2, sum))
mat2 <- rbind(sr, mat2)
mat2 <- mat2[ ,sort.list(mat2[1, ])]
mat2 <- mat2[2:(indnr+1), ]
mat1.sort <- matrix(NA, indnr, Ka)

for (i in 0:(Ka-2)) {
	difference <- matrix(NA, indnr, (Ka-i))
	for (j in 1:(Ka-i)) {
		difference[ ,j] <- abs(mat2[ ,(Ka-i)]-mat1[ ,j])
		}
	dsr <- rank(apply(difference, 2, sum))
	mat1 <- rbind(dsr, mat1)
	mat1 <- mat1[ ,sort.list(mat1[1, ])]
	mat1.sort[ ,(Ka-i)] <- mat1[2:(indnr+1), 1]
	mat1 <- mat1[2:(indnr+1), 2:(Ka-i)]
 	rm(difference, dsr)
	}
  mat1.sort[ ,1] <- mat1
  res <- array(data=NA, c(indnr, Ka, 2))
  res[ , , 1] <- mat1.sort
  res[ , , 2] <- mat2
  return(res)
}

# function which calculates coefficients of similarity for a series of structure runs with the same K
# after Nordborg et al. 2005 (symmetric similarity coefficients)

simil.structure <- function(rundata) {

d <- dim(rundata)
runnb <- d[3]
indnb <- d[1]
Ko <- d[2]

simil <- matrix(NA, runnb*(runnb-1)/2, 3)
no <- matrix(1/Ko, indnb, Ko)
s <- 1
frob <- function(x) {
	return(sqrt(sum(x^2))) }

for (m in 1:(runnb-1)) {
	for (n in (m+1):(runnb)) {
		pair <- ordermat(rundata[ , , m], rundata[ , , n], indnb, Ko)
		simil[s, 1] <- m
		simil[s, 2] <- n
		simil[s, 3] <- 1 - (frob(pair[ , , 1] - pair[ , , 2])/ sqrt( (frob(pair[ , , 1] - no)) * (frob(pair[ , , 2] - no))))
		s <- s+1
		}
	}

nam <- paste("simil", Ko, ".txt", sep="", collapse=NULL)
write.table(simil, file=nam, quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)		
return(simil)
}


# Main function to calculate and plot similarity coefficients
# After Nordborg et al. 2005

Structure.simil <- function(infile, pop=0, locprior=F, nbloc=pop) {

runtab <- read.table(infile, sep="\t", header=FALSE, row.names=NULL)
runtab <- runtab[sort.list(runtab[ ,1]), ]
runnumbers <- table(runtab[ ,1])
Knb <- length (runnumbers)

simils <- matrix(NA, Knb, 4)
z <- 1

if(scan(file=as.character(runtab[1, 2]), skip= 13, n=1, what="character")== "Run") hd <- 13
if(scan(file=as.character(runtab[1, 2]), skip= 14, n=1, what="character")== "Run") hd <- 14
if(scan(file=as.character(runtab[1, 2]), skip= 15, n=1, what="character")== "Run") hd <- 15

ran <- 0
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(file=as.character(runtab[1, 2]), skip= (hd+7), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(file=as.character(runtab[1, 2]), skip= (hd+8), n=1, what="character")== "RANDOMIZE") ran <- 1

startl <- (hd+23+ran)
modadm <- T
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "NO") {
  startl <- (hd+24+ran) 
  modadm <- F
  }
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "RECESSIVE") {
if(scan(file=as.character(runtab[1, 2]), skip= (hd+7), n=1, what="character")== "NO") {
  startl <- (hd+25+ran)
  modadm <- F } else startl <- (hd+24+ran)}

if (pop == 0) {
      
for (i in 1:Knb) {
	runfiles <- as.vector(runtab[z:(z+runnumbers[i]-1), 2])
	K <- scan(runfiles[1], skip=(hd+3), n=1)
	indnumb <- scan(runfiles[1], skip=(hd+1), n=1)

 	if (K != 1) {	
    	  modcor <- T
        if (modadm == T) { 
          if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 7), n=1, what="character")== "Allele") modcor <- F }
        if (modadm == F) { 
          if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 6), n=1, what="character")== "Allele") modcor <- F }

        startlines <- 0
        if (modadm == T & modcor == T) startlines <- startl + 3*K + 12
        if (modadm == T & modcor == F) startlines <- startl + 2*K + 12  #changed
        if (modadm == F & modcor == T) startlines <- startl + 3*K + 11
        if (modadm == F & modcor == F) startlines <- startl + 2*K + 11

  if (startlines == 0) print("The model you used to run Structure is not compatible with Structure2.2-sum")

		rundat <- array(data=NA, c(indnumb, K, runnumbers[i]))

		for (h in 1:runnumbers[i]) {
			rundat[ , , h] <- as.matrix(read.table(file=runfiles[h], nrows=indnumb, skip=startlines)[ ,5:(K+4)])
			}
		mode(rundat) <- "numeric"
		
		similarity <- simil.structure(rundat)
		simils[i, 1] <- K
		simils[i, 2] <- runnumbers[i]
		simils[i, 3] <- mean(similarity[ , 3])
		simils[i, 4] <- sd(similarity[ , 3]) }

	z <- z + runnumbers[i]
	rm(runfiles)
	}	}

if (locprior == T) startl <- startl + 1

if (pop > 0) {              
for (i in 1:Knb) {
		
	runfiles <- as.vector(runtab[z:(z+runnumbers[i]-1), 2])
	K <- scan(runfiles[1], skip=(hd+3), n=1)
	indnumb <- scan(runfiles[1], skip=(hd+1), n=1)

	if (K != 1) {	
        modcor <- T
        if (modadm == T) { 
          if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 6 + pop), n=1, what="character")== "Allele") modcor <- F }
        if (modadm == F) { 
          if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 5 + pop), n=1, what="character")== "Allele") modcor <- F }
        if(locprior == T) {   
          if (modadm == T) { 
            if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 6 + pop + nbloc + 2), n=1, what="character")== "Allele") modcor <- F }
          if (modadm == F) { 
           if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 5 + pop + nbloc + 2), n=1, what="character")== "Allele") modcor <- F }
          }
          
        startlines <- 0
        if(locprior == F) {
          if (modadm == T & modcor == T) startlines <- startl + 3*K + 11 + pop
          if (modadm == T & modcor == F) startlines <- startl + 2*K + 11 + pop
          if (modadm == F & modcor == T) startlines <- startl + 3*K + 10 + pop
          if (modadm == F & modcor == F) startlines <- startl + 2*K + 10 + pop
        }
        if(locprior == T) {
          if (modadm == T & modcor == T) startlines <- startl + 3*K + 11 + pop + nbloc + 3
          if (modadm == T & modcor == F) startlines <- startl + 2*K + 11 + pop + nbloc + 3
          if (modadm == F & modcor == T) startlines <- startl + 3*K + 10 + pop + nbloc + 3
          if (modadm == F & modcor == F) startlines <- startl + 2*K + 10 + pop + nbloc + 3
        }
       
  if (startlines == 0) print("The model you used to run structure is not compatible with Structure sum")

		rundat <- array(data=NA, c(indnumb, K, runnumbers[i]))

		for (h in 1:runnumbers[i]) {
			rundat[ , , h] <- as.matrix(read.table(file=runfiles[h], nrows=indnumb, skip=startlines)[ , 6:(K+5)])
			}
		mode(rundat) <- "numeric"
		
		similarity <- simil.structure(rundat)
		simils[i, 1] <- K
		simils[i, 2] <- runnumbers[i]
		simils[i, 3] <- mean(similarity[ , 3])
		simils[i, 4] <- sd(similarity[ , 3]) }

	z <- z + runnumbers[i]
	rm(runfiles)
	}	}
	
plot(simils[ ,1], (simils[ ,3]+ simils[ ,4]), type="p", pch=24, xlab="K", ylab="similarity coefficient", ylim=c(0, 1.25), axes=FALSE, frame=TRUE)
lines(simils[ ,1], simils[ ,3], type="b")
points(simils[ ,1], (simils[ ,3]- simils[ ,4]), pch=25)
axis(1, at=c(1:K))
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2))
simils <- rbind(c("K", "nb of runs", "mean similarity coeff.", "standard dev."), simils)
write.table(simils, file="simil_coefficients.txt", quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)		
}

### Function which plots the four diagnostic plots suggested by Evano et al. 2005

Structure.deltaK <- function(infile, pop=0, locprior=F) {
# OK 21.01.09
runtab <- read.table(infile, sep="\t", header=FALSE, row.names=NULL)
runtab <- runtab[sort.list(runtab[ ,1]), ]
runnumbers <- table(runtab[ ,1])
Knb <- length (runnumbers)

averages <- matrix(NA, Knb, 7)
z <- 1

if(scan(file=as.character(runtab[1, 2]), skip= 13, n=1, what="character")== "Run") hd <- 13
if(scan(file=as.character(runtab[1, 2]), skip= 14, n=1, what="character")== "Run") hd <- 14
if(scan(file=as.character(runtab[1, 2]), skip= 15, n=1, what="character")== "Run") hd <- 15

ran <- 0
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(file=as.character(runtab[1, 2]), skip= (hd+7), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(file=as.character(runtab[1, 2]), skip= (hd+8), n=1, what="character")== "RANDOMIZE") ran <- 1

startl <- (hd+23+ran)
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "NO") startl <- (hd+24+ran) 
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "RECESSIVE") {
if(scan(file=as.character(runtab[1, 2]), skip= (hd+7), n=1, what="character")== "NO") startl <- (hd+25+ran) else startl <- (hd+24+ran)}

if (pop==0) {
for (i in 1:Knb) {
	
	runfiles <- as.vector(runtab[z:(z+runnumbers[i]-1), 2])
	K <- scan(runfiles[1], skip=(hd+3), n=1)
	lnprobab <- matrix(NA, runnumbers[i], 2)

	for (r in 1: runnumbers[i]) {
		lnprobab[r, 1] <- K
		lnprobab[r, 2] <- as.numeric(scan(runfiles[r], skip=(startl + 3 +2*K), n=7, what="character")[7])
		}
		
	averages[i, 1] <- K
	averages[i, 2] <- runnumbers[i]
	averages[i, 3] <- mean(lnprobab[ , 2])
	averages[i, 4] <- sd(lnprobab[ , 2])
	rm(lnprobab)

	z <- z + runnumbers[i]
	rm(runfiles)
	}	 }

if (locprior == T) startl <- startl + 1

if (pop > 0) {         
for (i in 1:Knb) {
	
	runfiles <- as.vector(runtab[z:(z+runnumbers[i]-1), 2])
	K <- scan(runfiles[1], skip=(hd+3), n=1)
	lnprobab <- matrix(NA, runnumbers[i], 2)

	for (r in 1: runnumbers[i]) {
		lnprobab[r, 1] <- K
		lnprobab[r, 2] <- as.numeric(scan(runfiles[r], skip=(startl + 2 + 2*K + pop), n=7, what="character")[7])
		}
		
	averages[i, 1] <- K
	averages[i, 2] <- runnumbers[i]
	averages[i, 3] <- mean(lnprobab[ , 2])
	averages[i, 4] <- sd(lnprobab[ , 2])
	rm(lnprobab)

	z <- z + runnumbers[i]
	rm(runfiles)
	}	 }
	
for (i in 2:Knb) {
	averages[i, 5]<- averages[i, 3] - averages[(i-1), 3]
	}

for (i in 2:(Knb-1)) {
	averages[i, 6]<- abs(averages[(i+1), 5] - averages[i, 5])
	}

for (i in 2:(Knb-1)) {
	averages[i, 7]<- averages[i, 6] / averages[i, 4]
	}

averages.t <- rbind(c("K", "nb of runs", "mean LnPData", "standard dev.", "L'(K)", "L''(K)", "DeltaK"), averages) 
write.table(averages.t, file="DeltaK.txt", quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)		

par(mfrow=c(2, 2)) 
plot(averages[ ,1], averages[ ,3], type="p", xlab="K", ylab="Mean L(K)")
#error.bar(averages[ ,1], averages[ ,3], averages[3]-averages[4], averages[3] + averages[4])
plot(averages[ ,1], averages[ ,5], type="p", xlab="K", ylab="Mean L'(K)")
plot(averages[ ,1], averages[ ,6], type="p", xlab="K", ylab="Mean L''(K)")
plot(averages[ ,1], averages[ ,7], type="b", xlab="K", ylab="Mean DeltaK")
}

### This function adds a column to structure outputs with a nb for the groups, 
### which the individual was assigned to

Structure.groups <- function(infile, pop=0, locprior = F, nbloc=pop) {
if(scan(infile, skip= 13, n=1, what="character")== "Run") hd <- 13
if(scan(infile, skip= 14, n=1, what="character")== "Run") hd <- 14
if(scan(infile, skip= 15, n=1, what="character")== "Run") hd <- 15

ran <- 0
if(scan(infile, skip= (hd+6), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(infile, skip= (hd+7), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(infile, skip= (hd+8), n=1, what="character")== "RANDOMIZE") ran <- 1

startl <- (hd+23+ran)
modadm <- T
if(scan(infile, skip= (hd+6), n=1, what="character")== "NO") {
  startl <- (hd+24+ran) 
  modadm <- F
  }
if(scan(infile, skip= (hd+6), n=1, what="character")== "RECESSIVE") {
if(scan(infile, skip= (hd+7), n=1, what="character")== "NO") {
  startl <- (hd+25+ran)
  modadm <- F } else startl <- (hd+24+ran)}

K <- scan(infile, skip=(hd+3), n=1)
indnumb <- scan(infile, skip=(hd+1), n=1)
  
 if (pop == 0) {
    	  modcor <- T
        if (modadm == T) { 
          if(scan(infile, skip= (startl + 2*K + 7), n=1, what="character")== "Allele") modcor <- F }
        if (modadm == F) { 
          if(scan(infile, skip= (startl + 2*K + 6), n=1, what="character")== "Allele") modcor <- F }
        startlines <- 0
        if (modadm == T & modcor == T) startlines <- startl + 3*K + 12
        if (modadm == T & modcor == F) startlines <- startl + 2*K + 12
        if (modadm == F & modcor == T) startlines <- startl + 3*K + 11
        if (modadm == F & modcor == F) startlines <- startl + 2*K + 11

  if (startlines == 0) print("The model you used to run structure is not compatible with Structure sum")

	rundata <- as.matrix(read.table(file=infile, nrows=indnumb, skip=startlines)[ ,5:(K+4)])
	groups <- vector(mode="integer", indnumb)
 }

if (locprior == T) startl <- startl + 1
 
 if (pop > 0) {
        modcor <- T
        if (modadm == T) { 
          if(scan(infile, skip= (startl + 2*K + 6 + pop), n=1, what="character")== "Allele") modcor <- F }
        if (modadm == F) { 
          if(scan(infile, skip= (startl + 2*K + 5 + pop), n=1, what="character")== "Allele") modcor <- F }
        if(locprior == T) {   
          if (modadm == T) { 
            if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 6 + pop + nbloc + 2), n=1, what="character")== "Allele") modcor <- F }
          if (modadm == F) { 
           if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 5 + pop + nbloc + 2), n=1, what="character")== "Allele") modcor <- F }
          }

        startlines <- 0
        if (locprior == F){
          if (modadm == T & modcor == T) startlines <- startl + 3*K + 11 + pop
          if (modadm == T & modcor == F) startlines <- startl + 2*K + 11 + pop
          if (modadm == F & modcor == T) startlines <- startl + 3*K + 10 + pop
          if (modadm == F & modcor == F) startlines <- startl + 2*K + 10 + pop
        }
        if(locprior == T) {
          if (modadm == T & modcor == T) startlines <- startl + 3*K + 11 + pop + nbloc + 3
          if (modadm == T & modcor == F) startlines <- startl + 2*K + 11 + pop + nbloc + 3
          if (modadm == F & modcor == T) startlines <- startl + 3*K + 10 + pop + nbloc + 3
          if (modadm == F & modcor == F) startlines <- startl + 2*K + 10 + pop + nbloc + 3
        }


  if (startlines == 0) print("The model you used to run structure is not compatible with Structure sum")

	rundata <- as.matrix(read.table(file=infile, nrows=indnumb, skip=startlines)[ ,6:(K+5)])
	groups <- vector(mode="integer", indnumb)
 }
 
	for (i in 1:indnumb) {
		for (j in 1:K) {
			if (rundata[i, j]== max(rundata[i, ])) groups[i] <- j
		}}
	
	rundata <- cbind(rundata, groups)
	write.table(rundata, file="struct-groups.txt", quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)	
}	



## Function which orders the matrices resulting from a batch of Structure runs
## in a way that column 1 corresponds to the same cluster for each run
## The input for this function is a list of runs which were carried out for the same K.

Structure.order <- function(infile, pop=0, locprior = F, nbloc=pop) {
runtab <- read.table(infile, sep="\t", header=FALSE, row.names=NULL)
runtab <- runtab[sort.list(runtab[ ,1]), ]
runnumbers <- dim(runtab)[1]

if(scan(file=as.character(runtab[1, 2]), skip= 13, n=1, what="character")== "Run") hd <- 13
if(scan(file=as.character(runtab[1, 2]), skip= 14, n=1, what="character")== "Run") hd <- 14
if(scan(file=as.character(runtab[1, 2]), skip= 15, n=1, what="character")== "Run") hd <- 15

ran <- 0
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(file=as.character(runtab[1, 2]), skip= (hd+7), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(file=as.character(runtab[1, 2]), skip= (hd+8), n=1, what="character")== "RANDOMIZE") ran <- 1

startl <- (hd+23+ran)
modadm <- T
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "NO") {
  startl <- (hd+24+ran) 
  modadm <- F
  }
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "RECESSIVE") {
if(scan(file=as.character(runtab[1, 2]), skip= (hd+7), n=1, what="character")== "NO") {
  startl <- (hd+25+ran)
  modadm <- F } else startl <- (hd+24+ran)}

if (pop == 0) {
		
	runfiles <- as.vector(runtab[ , 2])
	K <- scan(runfiles[1], skip=(hd+3), n=1)
	indnumb <- scan(runfiles[1], skip=(hd+1), n=1)
	
    	  modcor <- T
        if (modadm == T) { 
          if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 7), n=1, what="character")== "Allele") modcor <- F }
        if (modadm == F) { 
          if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 6), n=1, what="character")== "Allele") modcor <- F }

        startlines <- 0
        if (modadm == T & modcor == T) startlines <- startl + 3*K + 12
        if (modadm == T & modcor == F) startlines <- startl + 2*K + 12
        if (modadm == F & modcor == T) startlines <- startl + 3*K + 11
        if (modadm == F & modcor == F) startlines <- startl + 2*K + 11

  if (startlines == 0) print("The model you used to run structure is not compatible with Structure sum")

		rundat <- array(data=NA, c(indnumb, K, runnumbers))
    
		for (h in 1:runnumbers) {
			rundat[ , , h] <- as.matrix(read.table(file=runfiles[h], nrows=indnumb, skip=startlines)[ ,5:(K+4)])
			}
    mode(rundat) <- "numeric"
    
		orderdat <- array(data=NA, c(indnumb, K, runnumbers))
		orderdat [ , , 1:2] <- ordermat(rundat[ , , 1], rundat[ , , 2], indnumb, K)
		if (runnumbers > 2) {
		  for (i in 3:runnumbers)
		  orderdat [ , , i] <- ordermat(rundat[ , , 1], rundat[ , , i], indnumb, K)[ , ,2]
			}
		}
		
if (locprior == T) startl <- startl + 1

if (pop > 0) {
		
	runfiles <- as.vector(runtab[ , 2])
	K <- scan(runfiles[1], skip=(hd+3), n=1)
	indnumb <- scan(runfiles[1], skip=(hd+1), n=1)

        modcor <- T
        if (modadm == T) { 
          if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 6 + pop), n=1, what="character")== "Allele") modcor <- F }
        if (modadm == F) { 
          if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 5 + pop), n=1, what="character")== "Allele") modcor <- F }
        if(locprior == T) {   
          if (modadm == T) { 
            if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 6 + pop + nbloc + 2), n=1, what="character")== "Allele") modcor <- F }
          if (modadm == F) { 
           if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 5 + pop + nbloc + 2), n=1, what="character")== "Allele") modcor <- F }
          }

        startlines <- 0     
        if (locprior == F){
          if (modadm == T & modcor == T) startlines <- startl + 3*K + 11 + pop
          if (modadm == T & modcor == F) startlines <- startl + 2*K + 11 + pop
          if (modadm == F & modcor == T) startlines <- startl + 3*K + 10 + pop
          if (modadm == F & modcor == F) startlines <- startl + 2*K + 10 + pop
        }
        if (locprior == T) {
          if (modadm == T & modcor == T) startlines <- startl + 3*K + 11 + pop + nbloc + 3
          if (modadm == T & modcor == F) startlines <- startl + 2*K + 11 + pop + nbloc + 3
          if (modadm == F & modcor == T) startlines <- startl + 3*K + 10 + pop + nbloc + 3
          if (modadm == F & modcor == F) startlines <- startl + 2*K + 10 + pop + nbloc + 3
        }

  if (startlines == 0) print("The model you used to run structure is not compatible with Structure sum")

		rundat <- array(data=NA, c(indnumb, K, runnumbers))

		for (h in 1:runnumbers) {
			rundat[ , , h] <- as.matrix(read.table(file=runfiles[h], nrows=indnumb, skip=startlines)[ , 6:(K+5)])
			}
		mode(rundat) <- "numeric"
		
		orderdat <- array(data=NA, c(indnumb, K, runnumbers))
		orderdat [ , , 1:2] <- ordermat(rundat[ , , 1], rundat[ , , 2], indnumb, K)
		if (runnumbers > 2) {
		  for (i in 3:runnumbers)
		  orderdat [ , , i] <- ordermat(rundat[ , , 1], rundat[ , , i], indnumb, K)[ , ,2]
			}
	}

for (i in 1:runnumbers) {
  nam <- paste(runfiles[i],"-ord.txt", sep="", collapse=NULL)
  write.table(orderdat[ , , i], file=nam, quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
}

averages <- apply(orderdat, c(1,2), mean)
sds <- apply(orderdat, c(1,2), sd)
maxmin <- apply(orderdat, c(1,2), max) - apply(orderdat, c(1,2), min)
write.table(averages, file="averages.txt", quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
write.table(sds, file="standard deviations.txt", quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
write.table(maxmin, file="max-min.txt", quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
}

## This function calculates a coefficient of clusteredness for each Structure run 
# according to Rosenberg et al. 2005. 

Structure.cluster <- function(infile, pop=0, locprior = F, nbloc=pop) {

runtab <- read.table(infile, sep="\t", header=FALSE, row.names=NULL)
runtab <- runtab[sort.list(runtab[ ,1]), ]
runnumbers <- dim(runtab)[1]

clusterness <- matrix(NA, runnumbers, 3)
clusterness[ , 2] <- runtab[ , 1]

if(scan(file=as.character(runtab[1, 2]), skip= 13, n=1, what="character")== "Run") hd <- 13
if(scan(file=as.character(runtab[1, 2]), skip= 14, n=1, what="character")== "Run") hd <- 14
if(scan(file=as.character(runtab[1, 2]), skip= 15, n=1, what="character")== "Run") hd <- 15

ran <- 0
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(file=as.character(runtab[1, 2]), skip= (hd+7), n=1, what="character")== "RANDOMIZE") ran <- 1
if(scan(file=as.character(runtab[1, 2]), skip= (hd+8), n=1, what="character")== "RANDOMIZE") ran <- 1

startl <- (hd+23+ran)
modadm <- T
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "NO") {
  startl <- (hd+24+ran) 
  modadm <- F
  }
if(scan(file=as.character(runtab[1, 2]), skip= (hd+6), n=1, what="character")== "RECESSIVE") {
if(scan(file=as.character(runtab[1, 2]), skip= (hd+7), n=1, what="character")== "NO") {
  startl <- (hd+25+ran)
  modadm <- F } else startl <- (hd+24+ran)}

if (pop == 0) {
		
	runfiles <- as.vector(runtab[ , 2])
	indnumb <- scan(runfiles[1], skip=(hd+1), n=1)
  
  for (i in 1:runnumbers) {            
  K <- scan(runfiles[i], skip=(hd+3), n=1)
  
   if (K != 1) {
        modcor <- T
        if (modadm == T) { 
          if(scan(file=as.character(runfiles[i]), skip= (startl + 2*K + 7), n=1, what="character")== "Allele") modcor <- F }
        if (modadm == F) { 
          if(scan(file=as.character(runfiles[i]), skip= (startl + 2*K + 6), n=1, what="character")== "Allele") modcor <- F }

        startlines <- 0
        if (modadm == T & modcor == T) startlines <- startl + 3*K + 12
        if (modadm == T & modcor == F) startlines <- startl + 2*K + 12
        if (modadm == F & modcor == T) startlines <- startl + 3*K + 11
        if (modadm == F & modcor == F) startlines <- startl + 2*K + 11

      rundat <- as.matrix(read.table(file=runfiles[i], nrows=indnumb, skip=startlines)[ ,5:(K+4)])
      mode(rundat) <- "numeric"
      
	   	divK <- matrix(1/K, indnumb, K)
		  clust1 <- (K/(K-1)) * apply((rundat - divK)^2, 1, sum)
		  clusterness [i, 3] <- sum(clust1) / indnumb
		} 
    if (K==1) clusterness [i, 3] <- "NA"
    clusterness [i, 1] <- runfiles[i]
  }	}

if (locprior == T) startl <- startl + 1

if (pop > 0) {
		
	runfiles <- as.vector(runtab[ , 2])
	indnumb <- scan(runfiles[1], skip=(hd+1), n=1)

  for (i in 1:runnumbers) {
  K <- scan(runfiles[i], skip=(hd+3), n=1)

	if (K!=1) {
        modcor <- T
        if (modadm == T) { 
          if(scan(file=as.character(runfiles[i]), skip= (startl + 2*K + 6 + pop), n=1, what="character")== "Allele") modcor <- F }
        if (modadm == F) { 
          if(scan(file=as.character(runfiles[i]), skip= (startl + 2*K + 5 + pop), n=1, what="character")== "Allele") modcor <- F }
        if(locprior == T) {   
          if (modadm == T) { 
            if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 6 + pop + nbloc + 2), n=1, what="character")== "Allele") modcor <- F }
          if (modadm == F) { 
           if(scan(file=as.character(runfiles[1]), skip= (startl + 2*K + 5 + pop + nbloc + 2), n=1, what="character")== "Allele") modcor <- F }
          }

        startlines <- 0
        if (locprior == F) {
         if (modadm == T & modcor == T) startlines <- startl + 3*K + 11 + pop
         if (modadm == T & modcor == F) startlines <- startl + 2*K + 11 + pop
         if (modadm == F & modcor == T) startlines <- startl + 3*K + 10 + pop
         if (modadm == F & modcor == F) startlines <- startl + 2*K + 10 + pop
        }
        if(locprior == T) {
          if (modadm == T & modcor == T) startlines <- startl + 3*K + 11 + pop + nbloc + 3
          if (modadm == T & modcor == F) startlines <- startl + 2*K + 11 + pop + nbloc + 3
          if (modadm == F & modcor == T) startlines <- startl + 3*K + 10 + pop + nbloc + 3
          if (modadm == F & modcor == F) startlines <- startl + 2*K + 10 + pop + nbloc + 3
        }

		rundat <- as.matrix(read.table(file=runfiles[i], nrows=indnumb, skip=startlines)[ , 6:(K+5)])
		mode(rundat) <- "numeric"
		
		divK <- matrix(1/K, indnumb, K)
		clust1 <- (K/(K-1)) * apply((rundat - divK)^2, 1, sum)
		clusterness [i, 3] <- sum(clust1) / indnumb
		}
   if (K==1) clusterness [i, 3] <- "NA" 
   clusterness [i, 1] <- runfiles[i] 
   }
	}
 
plot(clusterness[ ,2], clusterness[ ,3], xlab="K", ylab="Clusteredness", ylim=c(0,1), axes=FALSE, cex=1.8, frame=TRUE)
axis(1, at=c(1:K))
axis(2)

clusterness <- rbind(c("run", "K", "coef. of clusteredness"), clusterness)
write.table(clusterness, file="clustered_coefficients.txt", quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)		
}

