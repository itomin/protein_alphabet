
cle.dir.raw <- file.path(getwd(), "dev", "data", "pdbcle")


cle.files <- sort(list.files(cle.dir.raw))

seq.pdbs <- list()
seq.names <- list()


for(f in  cle.files){
  lines <- readLines(file.path(cle.dir.raw, f))
  seq <- lines[1]
  seq.pdbs <- seq
  pdb.name <- substr(f, 1,4)
  pdb.chain <- substr(f, 5,8)
  seq.name <- paste0(pdb.name, ":", pdb.chain)
  seqinr::write.fasta(sequences=seq.pdbs, names=seq.name, file.out=paste0(getwd(),"/dev/data/cath/", paste0(f, ".fasta")), as.string=T)
}




pairs.dir <- file.path(getwd(), "..", "data", "cath_interfold_pairs.csv")
df <- read.csv2(pairs.dir,  sep=";", header = F, col.names = c("tag"))
df

df.new <- df %>%
  separate(tag, c("C", "A", "T", "H", "PDB_ID", "PDB_CHAIN"), "\\.") %>% 
  select(PDB_ID, PDB_CHAIN)                


head(df.new)


write.table(df.new, file.path(getwd(),  "cath_interfold_pairs.csv"), quote = FALSE, row.names = FALSE, sep=";")
