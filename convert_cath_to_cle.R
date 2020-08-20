library(bio3d)
library(readtext)
library(dplyr)


pdb.dir <- file.path(getwd(), "..", "data", "dompdb")
cle.dir.raw <- file.path(getwd(), "..", "data", "pdbcle")

read.pdb(file = "12as")

pdb.files <- sort(list.files(pdb.dir))

for(f in  pdb.files){
  pdb.file <- file.path(pdb.dir, f)
  cle.file <- file.path(cle.dir.raw, f)
  system(paste("./pdbcle", pdb.file, cle.file))
  
  cle.raw <- readtext(cle.file)
  cle.raw <- iconv(cle.raw$text, "UTF-8", "UTF-8",sub='')
  cle.raw <- gsub("[^[:alnum:][:space:]]","",cle.raw)
  cle.raw <- strsplit(cle.raw, "\n")[[1]]
  cle.raw <- cle.raw[(length(cle.raw) - 1) : length(cle.raw)]
  write(cle.raw, paste0(cle.file))
}
