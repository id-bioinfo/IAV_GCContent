using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SummarySVMTrain
{
    class SummarySVMTrain
    {
        //select the best model based on the cross-validation performance
        //and save the best performance to file, svm_bestmodel_XXXX.csv
        static void Main(string[] args)
        {
            string prefix_path = "/Scripts/SVM/12_11To1_lineage/";
            List<List<List<double>>> svm_bestmodel_repeat_str = new List<List<List<double>>>();
            List<List<List<double>>> svm_bestmodel_crossvalidte = new List<List<List<double>>>();
            for (int i = 0; i < 12; ++i)
            {
                string accuracylog_file = prefix_path + "result_12_11To1_lineage_c_5fold_" + i + ".txt";
                List<String[]> ls = ReadTSV(accuracylog_file);
                for (int j = 0; j < ls.Count; ++j)
                {
                    if (ls[j][0].Contains("train_test"))
                    {
                        svm_bestmodel_repeat_str.Add(new List<List<double>>());
                    }
                    else
                    {
                        List<double> temp = new List<double>();
                        for (int k = 0; k < ls[j].Length; ++k)
                            if (!String.IsNullOrWhiteSpace(ls[j][k]) && !String.IsNullOrEmpty(ls[j][k]))
                                temp.Add(Double.Parse(ls[j][k]));
                        svm_bestmodel_repeat_str.Last().Add(temp);
                    }
                }

                string validatelog_file = prefix_path + "validate_12_11To1_lineage_c_5fold_" + i + ".txt";
                List<string> validatelog = ReadLines(validatelog_file);
                List<List<double>> validateles = new List<List<double>>();
                for (int j = 0; j < validatelog.Count; ++j)
                {
                    if(validatelog[j].Contains("Cross Validation Accuracy"))
                    {
                        validateles.Add(new List<double>());
                        validateles.Last().Add(Double.Parse(validatelog[j].Replace("%","").Split(' ').Last()));
                    }
                    else if(validatelog[j].Contains("Cross Validation Balanced Accuracy"))
                    {
                        validateles.Last().Add(Double.Parse(validatelog[j].Replace("%", "").Split(' ').Last()));
                    }
                }
                svm_bestmodel_crossvalidte.Add(validateles);
                Console.WriteLine(svm_bestmodel_repeat_str.Last().Count + "," + validateles.Count);
            }


            List<List<double>> svm_bestmodel_repeats = new List<List<double>>();
            for (int i = 0; i < svm_bestmodel_repeat_str.Count; ++i)
            {
                double best_balance_accuracy = 0, balance_accuracy;
                int bestIdx = 0;
                double accuracy_train, positive_accuracy_train, negative_accuracy_train;
                double accuracy_test, positive_accuracy_test, negative_accuracy_test;
                double avian_accuracy_test, lineage_Sw2P_accuracy_test;
                double validate_balance_accuracy;
                double validate_accuracy;
                double nweight, pweight, kernal, cost;
 
                int myIdx = 0;
                for (int j = 0; j < svm_bestmodel_repeat_str[i].Count; j += 5)
                {
                    if (svm_bestmodel_repeat_str[i][j].Count == 10 && svm_bestmodel_repeat_str[i][j + 1].Count == 10)
                    {
                        myIdx = j;
                        positive_accuracy_train = svm_bestmodel_repeat_str[i][myIdx][2] / svm_bestmodel_repeat_str[i][myIdx][3];
                        negative_accuracy_train = svm_bestmodel_repeat_str[i][myIdx][4] / svm_bestmodel_repeat_str[i][myIdx][5];
                        ++myIdx;
                        positive_accuracy_test = svm_bestmodel_repeat_str[i][myIdx][2] / svm_bestmodel_repeat_str[i][myIdx][3];
                        negative_accuracy_test = svm_bestmodel_repeat_str[i][myIdx][4] / svm_bestmodel_repeat_str[i][myIdx][5];
                        ++myIdx;
                        avian_accuracy_test = svm_bestmodel_repeat_str[i][myIdx][0] / svm_bestmodel_repeat_str[i][myIdx][1];
                        ++myIdx;
                        lineage_Sw2P_accuracy_test = svm_bestmodel_repeat_str[i][myIdx][0] / svm_bestmodel_repeat_str[i][myIdx][1];
                        validate_accuracy = svm_bestmodel_crossvalidte[i][j / 5][0];
                        validate_balance_accuracy = svm_bestmodel_crossvalidte[i][j / 5][1];
                        balance_accuracy = validate_balance_accuracy; 
                        if (balance_accuracy > best_balance_accuracy)
                        {
                            bestIdx = j;
                            best_balance_accuracy = balance_accuracy;
                        }
                    }
                }
                List<double> one_repeat = new List<double>();
                myIdx = bestIdx;
                accuracy_train = svm_bestmodel_repeat_str[i][myIdx][0] / svm_bestmodel_repeat_str[i][myIdx][1];
                positive_accuracy_train = svm_bestmodel_repeat_str[i][myIdx][2] / svm_bestmodel_repeat_str[i][myIdx][3];
                negative_accuracy_train = svm_bestmodel_repeat_str[i][myIdx][4] / svm_bestmodel_repeat_str[i][myIdx][5];
                ++myIdx;
                accuracy_test = svm_bestmodel_repeat_str[i][myIdx][0] / svm_bestmodel_repeat_str[i][myIdx][1];
                positive_accuracy_test = svm_bestmodel_repeat_str[i][myIdx][2] / svm_bestmodel_repeat_str[i][myIdx][3];
                negative_accuracy_test = svm_bestmodel_repeat_str[i][myIdx][4] / svm_bestmodel_repeat_str[i][myIdx][5];
                ++myIdx;
                avian_accuracy_test = svm_bestmodel_repeat_str[i][myIdx][0] / svm_bestmodel_repeat_str[i][myIdx][1];
                ++myIdx;
                lineage_Sw2P_accuracy_test = svm_bestmodel_repeat_str[i][myIdx][0] / svm_bestmodel_repeat_str[i][myIdx][1];
                ++myIdx;
                lineage_Sw2N_accuracy_test = svm_bestmodel_repeat_str[i][myIdx][0] / svm_bestmodel_repeat_str[i][myIdx][1];

                nweight = svm_bestmodel_repeat_str[i][bestIdx][6];
                pweight = svm_bestmodel_repeat_str[i][bestIdx][7];
                kernal = svm_bestmodel_repeat_str[i][bestIdx][8];
                cost = svm_bestmodel_repeat_str[i][bestIdx][9];

                one_repeat.Add(accuracy_train);
                one_repeat.Add(positive_accuracy_train);
                one_repeat.Add(negative_accuracy_train);
                one_repeat.Add((positive_accuracy_train + negative_accuracy_train) / 2);
                one_repeat.Add(accuracy_test);
                one_repeat.Add(positive_accuracy_test);
                one_repeat.Add(negative_accuracy_test);
                one_repeat.Add((positive_accuracy_test + negative_accuracy_test) / 2);
                one_repeat.Add(avian_accuracy_test);
                one_repeat.Add(lineage_Sw2P_accuracy_test);
      
                one_repeat.Add(nweight);
                one_repeat.Add(pweight);
                one_repeat.Add(kernal);
                one_repeat.Add(cost);

                svm_bestmodel_repeats.Add(one_repeat);
            }

            string outputfile = prefix_path + "svm_bestmodel_12_11To1_lineage.csv";

            StreamWriter fileWriter = new StreamWriter(outputfile, false, Encoding.ASCII);
            fileWriter.WriteLine("index,accuracy_train,positive_accuracy_train,negative_accuracy_train,balance_accuracy_train,accuracy_test,positive_accuracy_test,negative_accuracy_test,balance_accuracy_test,avian_accuracy,earliest,nweight,pweight,kernal,c");
            for (int i = 0; i < svm_bestmodel_repeats.Count; ++i)
            {
                fileWriter.Write(i);
                for (int j = 0; j < svm_bestmodel_repeats[i].Count; ++j)
                    fileWriter.Write("," + svm_bestmodel_repeats[i][j]);
                fileWriter.Write("\n");
            }
            fileWriter.Close();
        }

        static List<String[]> ReadCSV(string filePathName)
        {
            List<String[]> ls = new List<String[]>();
            StreamReader fileReader = new StreamReader(filePathName, Encoding.Default);
            string line = "";
            while ((line = fileReader.ReadLine()) != null)
            {
                line = line.Trim();
                if (line.Length > 0)
                {
                    ls.Add(line.Split(','));
                }
            }
            fileReader.Close();
            return ls;
        }
        static List<String[]> ReadTSV(string filePathName)
        {
            List<String[]> ls = new List<String[]>();
            StreamReader fileReader = new StreamReader(filePathName, Encoding.Default);
            string line = "";
            while ((line = fileReader.ReadLine()) != null)
            {
                line = line.Trim();
                if (line.Length > 0)
                {
                    ls.Add(line.Split('\t'));
                }
            }
            fileReader.Close();
            return ls;
        }

        static List<string> ReadLines(string filePathName)
        {
            List<string> ls = new List<string>();
            StreamReader fileReader = new StreamReader(filePathName, Encoding.Default);
            string line = "";
            while ((line = fileReader.ReadLine()) != null)
            {
                line = line.Trim();
                if (line.Length > 0)
                {
                    ls.Add(line);
                }
            }
            fileReader.Close();
            return ls;
        }
    }
}
