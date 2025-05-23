#!/bin/sh

cd ~/yytao/influenza/data/filter
method=DOWNPASS
gene=PB2
intree=$gene"_nt_pal2nal.fasttree_midpoint_collapsed_InnodeNameAdded"
hostinfofile=$gene"_host_info.tsv"
cd pastml
outfolder=$gene"_"$method
mkdir $outfolder

pastml --tree ../$intree --data $hostinfofile --columns Host --html_compressed $outfolder/$gene"_host_map.html" --prediction_method $method --threads 10 --work_dir $outfolder --pajek $outfolder/$gene"_host_map.pajek" --tip_size_threshold 40