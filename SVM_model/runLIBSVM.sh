#!/bin/sh

#compute the features of GC content and GC dinucleotide frequency (CpG, GpC, GpG and CpC) from fasta files in LIBSVM format
python computeFeature.py
#apply our SVM model for the risk assessment
modelfile=cds/training_svm.model
testfile=demo_minkH5_cds/gc_content_dinucleotide_libsvm.txt
#the result is saved in $testfile"_predict"
#the installation and usage of libsvm-3.3 please refer to https://github.com/cjlin1/libsvm.
libsvm-3.3/svm-predict $testfile $modelfile $testfile"_predict"

#traingfile and testfile are in folder "SVM_model"
#training
libsvm-3.3/svm-train -t 0 -w-1 16 $trainfile $modelfile
#testing
libsvm-3.3/svm-predict $testfile $modelfile $testfile"_predict"