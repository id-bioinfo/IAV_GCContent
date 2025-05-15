library(pcalg)

##test gc against 4 dinucleotides for 12 mammalian lineages, Hu1-4, Sw1-5, Eq, Ca1-2
gc <- read.table('D:\\workplace\\projects\\Influenza evolution\\gc_content\\genome\\gc_Hu1_year.txt', header = F, sep = "\t")
gc <- gc[,2]
cpg <- read.table('D:\\workplace\\projects\\Influenza evolution\\dinucleotide\\CpG_Hu1_year.txt', header = F, sep = "\t")
cpg <- cpg[,2]
gpc <- read.table('D:\\workplace\\projects\\Influenza evolution\\dinucleotide\\GpC_Hu1_year.txt', header = F, sep = "\t")
gpc <- gpc[,2]
gpg <- read.table('D:\\workplace\\projects\\Influenza evolution\\dinucleotide\\GpG_Hu1_year.txt', header = F, sep = "\t")
gpg <- gpg[,2]
cpc <- read.table('D:\\workplace\\projects\\Influenza evolution\\dinucleotide\\CpC_Hu1_year.txt', header = F, sep = "\t")
cpc <- cpc[,2]

X <- cbind(gc,cpg)
res <- lingam(X)
cat("CpG estimated DAG:\n")
res$Bpruned

X <- cbind(gc,gpc)
res <- lingam(X)
cat("GpC estimated DAG:\n")
res$Bpruned

X <- cbind(gc,gpg)
res <- lingam(X)
cat("GpG estimated DAG:\n")
res$Bpruned

X <- cbind(gc,cpc)
res <- lingam(X)
cat("CpC estimated DAG:\n")
res$Bpruned

