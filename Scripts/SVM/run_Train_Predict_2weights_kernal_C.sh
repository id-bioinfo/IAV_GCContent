#!/bin/bash 
libsvm_path=/media/data/ytye/package/libsvm-3.3
svm_train=$libsvm_path/svm-train
svm_predict=$libsvm_path/svm-predict
predixPath=12_11To1_lineage/dataset_12_11To1Lineage

for((i=0;i<12;i++)); 
do
validate_log=validate_12_11To1_lineage_c_5fold_"$i".txt
result_log=result_12_11To1_lineage_c_5fold_"$i".txt
rm $validate_log
rm $result_log
trainfile=$predixPath/train_dataset_$i
modelfile=$predixPath/model_$i
echo "train_test_"$i >> $result_log
for((n_weight_idx=1;n_weight_idx<=5;n_weight_idx++)); 
do
    n_weight=$(echo "$n_weight_idx * 2 - 1" | bc)
    for((p_weight_idx=1;p_weight_idx<=5;p_weight_idx++)); 
    do
        p_weight=$(echo "$p_weight_idx * 2 - 1" | bc)
        for((kernal_t=0;kernal_t<=3;kernal_t++)); 
        do
            for((cost_idx=1;cost_idx<=6;cost_idx++)); 
            do
                if [ "$cost_idx" -eq 1 ]; then
                    cost=0.125
                else
                    cost=$(echo "$cost * 2" | bc)
                fi
                $svm_train -t $kernal_t -w-1 $n_weight -w1 $p_weight -c $cost -v 5 $trainfile >> $validate_log
                $svm_train -t $kernal_t -w-1 $n_weight -w1 $p_weight -c $cost $trainfile $modelfile >> $validate_log

                testfile=$trainfile
                $svm_predict $testfile $modelfile $testfile"_predict"
                awk '{print $1}' $testfile > $testfile"_class"
                paste -d "\t" $testfile"_class" $testfile"_predict" > $testfile"_class_predict"
                awk -F '\t' -v cost=$cost -v n_weight=$n_weight -v p_weight=$p_weight -v kernal_t=$kernal_t '{if($1==$2){tp++; if($1==1) p_tp++; if($1==-1) n_tp++;}; if($1==1) p_count++; if($1==-1) n_count++;} END{print tp"\t"NR"\t"p_tp"\t"p_count"\t"n_tp"\t"n_count"\t"n_weight"\t"p_weight"\t"kernal_t"\t"cost}' $testfile"_class_predict" >> $result_log

                testfile=$predixPath/test_dataset_$i
                $svm_predict $testfile $modelfile $testfile"_predict"
                awk '{print $1}' $testfile > $testfile"_class"
                paste -d "\t" $testfile"_class" $testfile"_predict" > $testfile"_class_predict"
                awk -F '\t' -v cost=$cost -v n_weight=$n_weight -v p_weight=$p_weight -v kernal_t=$kernal_t '{if($1==$2){tp++; if($1==1) p_tp++; if($1==-1) n_tp++;}; if($1==1) p_count++; if($1==-1) n_count++;} END{print tp"\t"NR"\t"p_tp"\t"p_count"\t"n_tp"\t"n_count"\t"n_weight"\t"p_weight"\t"kernal_t"\t"cost}' $testfile"_class_predict" >> $result_log
        
                testfile=$predixPath/av_dataset_$i
                $svm_predict $testfile $modelfile $testfile"_predict"
                awk '{print $1}' $testfile > $testfile"_class"
                paste -d "\t" $testfile"_class" $testfile"_predict" > $testfile"_class_predict"
                awk -F '\t' -v cost=$cost -v n_weight=$n_weight -v p_weight=$p_weight -v kernal_t=$kernal_t '{if($1==$2){tp++; if($1==1) p_tp++; if($1==-1) n_tp++;}; if($1==1) p_count++; if($1==-1) n_count++;} END{print tp"\t"NR"\t"p_tp"\t"p_count"\t"n_tp"\t"n_count"\t"n_weight"\t"p_weight"\t"kernal_t"\t"cost}' $testfile"_class_predict" >> $result_log   

                #testfile=$predixPath/earliest_dataset_$i
                #$svm_predict $testfile $modelfile $testfile"_predict"
                #awk '{print $1}' $testfile > $testfile"_class"
                #paste -d "\t" $testfile"_class" $testfile"_predict" > $testfile"_class_predict"
                #awk -F '\t' -v cost=$cost -v n_weight=$n_weight -v p_weight=$p_weight -v kernal_t=$kernal_t '{if($1==$2){tp++; if($1==1) p_tp++; if($1==-1) n_tp++;}; if($1==1) p_count++; if($1==-1) n_count++;} END{print tp"\t"NR"\t"p_tp"\t"p_count"\t"n_tp"\t"n_count"\t"n_weight"\t"p_weight"\t"kernal_t"\t"cost}' $testfile"_class_predict" >> $result_log
            done
        done
    done
done
done
