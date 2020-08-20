library(seqinr)
library(dplyr)
library(tidyr)
library(splitstackshape)

trans.matrix <- function(X, prob=T)
{
  tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
  if(prob) tt <- tt / rowSums(tt)
  tt
}


anno <- read.csv(file="../data/cath_annotation.csv", sep=";", header = F) %>% 
        separate(V2, c("C", "A", "TO", "H"), sep="\\.") %>% rename(seq = V1)


cle <- read.csv(file="../data/cath_pdb_cle_v001.csv", sep=",", header = F, stringsAsFactors = F) 
cle <- cle %>% mutate(L = nchar(V2))
sizes <- cle$L

cle.sep <- cSplit(cle, splitCols = "V2", sep="", stripWhite = F)
cle.sep <- cle.sep %>% select(-V1, -L)

trans.prob <- trans.matrix(as.matrix(cle.sep))
counts <- table(as.matrix(cle.sep))
initialprobs <- counts / sum(counts)
alphabet <- names(initialprobs)


random.mc <- function(alphabet, initialprobs, trans.prob, seqlength){
  new.sequence <- character()
  first.cle <- sample(alphabet, 1, rep=TRUE, prob=initialprobs)
  new.sequence[1] <- first.cle 
  for (i in 2:seqlength){
      prev.cle <- new.sequence[i-1]     # Get the previous nucleotide in the new sequence
      probabilities <- trans.prob[prev.cle,]
      next.cle <- sample(alphabet, 1, rep=TRUE, prob=probabilities)
      new.sequence[i] <- next.cle 
  }
  return(new.sequence)
}


for(i in 1:100){
  df.random <- do.call(rbind, lapply(1:nrow(cle), function(i){
    L <- random.mc(alphabet, initialprobs, trans.prob, sample(sizes, 1))
    paste(L, collapse = "")
  }))
  
  write.csv(df.random, paste0("random/random_mc_", i, ".csv"), row.names = F, col.names = F)
}
warnings()
