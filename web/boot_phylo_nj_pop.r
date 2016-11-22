## Bootstrap of NJ tree of populations. As input takes bootstraped tree, source data in geneind format and number of permutations.

boot.phylo.nj.pop <- function(origtree, genindata, nperm) {

# Custom variant of function sample to correct its behaviour when number of individuals within respective population is 1.
  mysample <- function(x){
    if (length(x) == 1) {out <- x}
      else {out <- sample(x, replace=TRUE)}
      return(out)
    }

# Construction of NJ tree from bootstraped dataset for further comparison with original tree.
## Here it is possible to adjust type and/or parameters of distance matrix and/or tree constructioning.
  treemaker <- function (ind, listgenpopobj) {
    OpravaMatice <- dist.genpop(listgenpopobj[[ind]], method=1, diag=TRUE, upper=TRUE)
    OpravaMatice[OpravaMatice == "Inf"] <- 10
    nj(OpravaMatice)
    }

# Extraction of needed information from orginal data and preparation of temporal variables.
  genindata2 <- genindata
  data <- genindata@tab
  pop <- genindata@pop
  indiv <- 1:nrow(genindata@tab)
  res <- vector(mode="list", length=nperm)

# In every step (from 1 to number of permutations) sample within respective populations of source data and write it into list of geneind objects.
  for (index in 1:nperm) {
    struct <- unlist(tapply(indiv, pop, mysample))
    new.tab <- data[struct,]
    genindata2@tab <- new.tab
# For correct writing of modified data to the list.
    res[[index]] <- genindata2
    }

# Conversion of all genind objects in the list to genpop objects.
  population <- lapply(res, genind2genpop)
# Construction of trees and saving of them into the list of trees.
  trees <- vector("list", nperm)
  trees <- lapply(1:nperm, treemaker, population)
# Compare every bootstraped tree with original one and return vector of bootstrap support values.
  for (i in 1:nperm) storage.mode(trees[[i]]$Nnode) <- "integer"
  storage.mode(origtree$Nnode) <- "integer"
  pp <- prop.part(trees)
  pp <- postprocess.prop.part(pp)
  ans <- prop.clades(origtree, part=pp, rooted=FALSE)
  return(ans)

  }
