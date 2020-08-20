
cle.dir.raw <- file.path(getwd(), "..", "data", "pdbcle")

cle.files <- sort(list.files(cle.dir.raw))

seq.cle <- list()
seq.names <- list()

for(f in  cle.files){
  lines <- readLines(file.path(cle.dir.raw, f))
  cle <- lines[2]
  seq.cle <- c(seq.cle, cle)
  pdb.name <- substr(f, 1,4)
  pdb.chain <- substr(f, 5,8)
  seq.names <- c(seq.names, paste0(pdb.name, ":", pdb.chain))
}

seqinr::write.fasta(sequences=seq.cle, names=seq.names, file.out=paste0(getwd(),"/","cath_pdb_cle_v001.fasta"), as.string=T)


pairs.dir <- file.path(getwd(), "..", "data", "cath_interfold_pairs.csv")
df <- read.csv2(pairs.dir,  sep=";", header = F, col.names = c("tag"))
df

df.new <- df %>%
separate(tag, c("C", "A", "T", "H", "PDB_ID", "PDB_CHAIN"), "\\.") %>%
select(PDB_ID, PDB_CHAIN)


head(df.new)


write.table(df.new, file.path(getwd(),  "cath_interfold_pairs.csv"), quote = FALSE, row.names = FALSE, sep=";")
