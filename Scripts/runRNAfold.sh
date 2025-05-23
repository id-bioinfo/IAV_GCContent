#!/bin/sh

cd ~/yytao/influenza/data/filter
gene=PB2
infile=$gene"_nt.pal2nal"
outfile=./rnafold/$gene"_rnafold.txt"

~/yytao/package/ViennaRNA-2.6.4/src/bin/RNAfold -j20 --MEA --noPS $infile > $outfile
