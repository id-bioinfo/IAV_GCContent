#!/bin/sh

cd ./filter

arraygene=(HA NA NP PA PB1 PB2 M1 NS1)
for gene in ${arraygene[@]}
do
  echo $gene
  aafile=$gene"_aa_msa_all.fasta"
  ntfile=../$gene"_nt.fasta"
  outfile=../$gene"_nt.pal2nal"
  ~/yytao/package/pal2nal/pal2nal.pl $aafile $ntfile -output fasta > $outfile
done
