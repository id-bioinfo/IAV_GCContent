#!/bin/sh

gene=PB2
cd ~/yytao/influenza/data/filter/partition/$gene

for((i=0;i<14;i++));
do
  infile=$gene"_aa_"$i".fasta"
  outfile=$gene"_aa_msa_"$i".fasta"
  ~/yytao/package/MAFFT/mafft --auto --thread 16 $infile > $outfile
done
