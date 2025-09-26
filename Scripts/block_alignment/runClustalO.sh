#!/bin/sh

gene=PB2
cd ./filter/partition/$gene

##the order of alignments is based on the guide tree from the consensus sequences of the 14 sub-MSAs
##orders are different for different genes

~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 $gene"_aa_msa_0.fasta" --p2 $gene"_aa_msa_1.fasta" -o output_h3_1
~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 $gene"_aa_msa_2.fasta" --p2 output_h3_1 -o output_h3

~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 $gene"_aa_msa_3.fasta" --p2 $gene"_aa_msa_4.fasta" -o output_h1_1
~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 output_h1_1 --p2 $gene"_aa_msa_5.fasta" -o output_h1

~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 $gene"_aa_msa_6.fasta" --p2 $gene"_aa_msa_10.fasta" -o output_1
~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 output_1 --p2 $gene"_aa_msa_13.fasta" -o output_2
~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 output_2 --p2 $gene"_aa_msa_9.fasta" -o output_3
~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 output_3 --p2 $gene"_aa_msa_7.fasta" -o output_4

~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 $gene"_aa_msa_11.fasta" --p2 $gene"_aa_msa_12.fasta" -o output_5
~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 output_4 --p2 output_5 -o output_6
~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 output_6 --p2 $gene"_aa_msa_8.fasta" -o output_7

~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 output_h1 --p2 output_7 -o output_8

~/yytao/package/clustalo/clustalo-1.2 --auto --threads=16 --p1 output_8 --p2 output_h3 -o ../../$gene"_aa_msa_all.fasta"
