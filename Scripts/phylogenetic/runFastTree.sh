#!/bin/sh

cd ./filter
gene=PB2
infile=$gene"_nt.pal2nal"
outfile=$gene"_nt_pal2nal.fasttree"

~/yytao/package/FastTree/FastTreeDouble -nt -gtr -gamma $infile > $outfile
