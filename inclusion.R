source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("cleaver")
library(cleaver)
library(Biostrings)

fasta <- readAAStringSet(file = "")

sequences <- as.character(fasta)

tryptic.peptides <- cleave(sequences, enzym = "trypsin")

are_peps_MSok <- function (x) {
  numAA <- nchar(x)
  y <- x[numAA >= 6 & numAA < 35]
  y <- y[grepl("R|K", y)]
  return(y)
}

msOK.peptides <- lapply(tryptic.peptides, are_peps_MSok)
unique.msOK.peptides <- unique(unlist(msOK.peptides))

aa_masses <- read.delim("amino acid masses (acc5)")
water.mass <- (1.00794*2)+15.9994

calc_mass <- function (this.seq) {
  indices <- match(unlist(strsplit(this.seq,"")), aa_masses$X1.letter.code)
  return(round(sum(aa_masses$Monoisotopic[indices]) + water.mass, 5))
}

unique.msOK.peptides.masses <- sapply(as.character(unique.msOK.peptides), calc_mass)

proton.mass <- 1.007276

singly.charged.masses <- unique.msOK.peptides.masses + proton.mass

doubly.charged.masses <- (singly.charged.masses + proton.mass)/2

all.masses <- c(singly.charged.masses, doubly.charged.masses)

df <- data.frame(mass = all.masses, formula = "", formula_type = "", species = "", polarity = "Positive", start = "", end = "", nce = "", nce.type = "", msx = "", Comment = names(all.masses))
df.export <- df
inclusion.template.colnames <- c("Mass [m/z]", "Formula [M]", "Formula type", "Species	CS [z]", "Polarity", "Start [min]", "End [min]", "(N)CE", "(N)CE type", "MSX ID", "Comment")

colnames(df.export) <- inclusion.template.colnames
write.csv(df.export, "inclusion_list.csv", quote = F, row.names = F)
