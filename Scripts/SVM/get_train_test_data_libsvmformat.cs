using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace change2libsvm_format
{
    class change2libsvm_format
    {
        static void Main_66_10To2_lineage(string[] args) //Main_66_10To2_lineage
        {
            string allfeature_file = "../genomic_features/all_features_seqname_lineage_year.tsv";
            List<String[]> ls = ReadTSV(allfeature_file);
            List<String[]> pe_ls = new List<String[]>();
            List<String[]> sp_ls = new List<String[]>();
            List<String[]> av_ls = new List<String[]>();
            HashSet<string> pe_lineages = new HashSet<string>(new string[] { "Hu1", "Hu2", "Hu3", "Hu4", "Eq", "Ca1", "Ca2", "Sw1", "Sw2", "Sw3", "Sw4", "Sw5" });
            for (int i = 0; i < ls.Count; ++i)
            {
                if (pe_lineages.Contains(ls[i][42])) pe_ls.Add(ls[i]);
                else if (ls[i][42] == "Av") av_ls.Add(ls[i]);
                else sp_ls.Add(ls[i]);
            }
            //pe_lineage2features
            Dictionary<string, List<String[]>> lineage_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < pe_ls.Count; ++i)
            {
                string lineage = pe_ls[i][42];
                if (!lineage_features.ContainsKey(lineage))
                {
                    lineage_features.Add(lineage, new List<String[]>());
                }
                lineage_features[lineage].Add(pe_ls[i]);
            }
            Console.WriteLine(lineage_features.Count);

            //list of lineage features
            List<List<String[]>> lineage_feature_list = new List<List<String[]>>();
            foreach (var item in lineage_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                lineage_feature_list.Add(sort_list);
            }
            Console.WriteLine(lineage_feature_list.Count);

            //sporadic
            //sp2features
            Dictionary<string, List<String[]>> sp_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < sp_ls.Count; ++i)
            {
                string sp = sp_ls[i][42];
                if (!sp_features.ContainsKey(sp))
                {
                    sp_features.Add(sp, new List<String[]>());
                }
                sp_features[sp].Add(sp_ls[i]);
            }
            Console.WriteLine(sp_features.Count);

            //list of sp features
            List<String[]> sp_feature_list = new List<String[]>();
            foreach (var item in sp_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                sp_feature_list.AddRange(sort_list);
            }
            Console.WriteLine(sp_feature_list.Count);

            //avian
            //av2features
            Dictionary<string, List<String[]>> av_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < av_ls.Count; ++i)
            {
                string av = av_ls[i][41].Split('|')[2]; //AvHxNx
                if (!av_features.ContainsKey(av))
                {
                    av_features.Add(av, new List<String[]>());
                }
                av_features[av].Add(av_ls[i]);
            }
            Console.WriteLine(av_features.Count);

            //list of av features
            List<String[]> av_feature_list = new List<String[]>();
            foreach (var item in av_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                av_feature_list.AddRange(sort_list);
            }
            Console.WriteLine(av_feature_list.Count);

            string earliest_lineage_file = "../genomic_features/earliest_features.tsv";
            List<String[]> earliest_feature_list = ReadTSV(earliest_lineage_file);

            int datasetId = 0;
            Random rnd = new Random();
            for (int lineage_index_1 = 0; lineage_index_1 < lineage_feature_list.Count; lineage_index_1++)
            {
                for (int lineage_index_2 = lineage_index_1 + 1; lineage_index_2 < lineage_feature_list.Count; lineage_index_2++)
                {
                    string trainfile = "66_10To2_lineage/dataset_66_10To2Lineage/train_dataset_" + datasetId;
                    string testfile = "66_10To2_lineage/dataset_66_10To2Lineage/test_dataset_" + datasetId;
                    string earliestfile = "66_10To2_lineage/dataset_66_10To2Lineage/earliest_dataset_" + datasetId;
                    string testfile_av = "66_10To2_lineage/dataset_66_10To2Lineage/av_dataset_" + datasetId;
                    StreamWriter trainfile_writer = new StreamWriter(trainfile, false, Encoding.ASCII);
                    StreamWriter testfile_writer = new StreamWriter(testfile, false, Encoding.ASCII);
                    StreamWriter earliestfile_writer = new StreamWriter(earliestfile, false, Encoding.ASCII);
                    StreamWriter testfile_av_writer = new StreamWriter(testfile_av, false, Encoding.ASCII);
                    datasetId++;
                    //positive
                    for (int j = 0; j < 12; ++j)
                    {
                        if (j != lineage_index_1 && j != lineage_index_2)
                        {
                            for (int m = 0; m < lineage_feature_list[j].Count; ++m)
                            {
                                trainfile_writer.Write("1");
                                for (int k = 1; k < 41; ++k)
                                {
                                    trainfile_writer.Write("\t" + k + ":" + lineage_feature_list[j][m][k]);
                                }
                                trainfile_writer.Write("\n");
                            }
                        }
                    }
                    for (int m = 0; m < lineage_feature_list[lineage_index_1].Count; ++m)
                    {
                        testfile_writer.Write("1");
                        for (int k = 1; k < 41; ++k)
                        {
                            testfile_writer.Write("\t" + k + ":" + lineage_feature_list[lineage_index_1][m][k]);
                        }
                        testfile_writer.Write("\n");
                    }
                    for (int m = 0; m < lineage_feature_list[lineage_index_2].Count; ++m)
                    {
                        testfile_writer.Write("1");
                        for (int k = 1; k < 41; ++k)
                        {
                            testfile_writer.Write("\t" + k + ":" + lineage_feature_list[lineage_index_2][m][k]);
                        }
                        testfile_writer.Write("\n");
                    }

                    //earliest sequences for test
                    earliestfile_writer.Write("1");
                    for (int k = 1; k < 41; ++k)
                    {
                        earliestfile_writer.Write("\t" + earliest_feature_list[lineage_index_1][k]);
                    }
                    earliestfile_writer.Write("\n");
                    earliestfile_writer.Write("1");
                    for (int k = 1; k < 41; ++k)
                    {
                        earliestfile_writer.Write("\t" + earliest_feature_list[lineage_index_2][k]);
                    }
                    earliestfile_writer.Write("\n");

    
                    //negative samples: 80% for training, 20% for testing
                    //sporadic
                    for (int m = 0; m < sp_feature_list.Count; ++m)
                    {
                        int sp_train_flag = rnd.Next(0, 10);  // creates a number between 0 and 9
                        if (sp_train_flag < 8)
                        {
                            trainfile_writer.Write("-1");
                            for (int k = 1; k < 41; ++k)
                            {
                                trainfile_writer.Write("\t" + k + ":" + sp_feature_list[m][k]);
                            }
                            trainfile_writer.Write("\n");
                        }
                        else
                        {
                            testfile_writer.Write("-1");
                            for (int k = 1; k < 41; ++k)
                            {
                                testfile_writer.Write("\t" + k + ":" + sp_feature_list[m][k]);
                            }
                            testfile_writer.Write("\n");
                        }
                    }
                    
                    //avian
                    for (int m = 0; m < av_feature_list.Count; ++m)
                    {
                        int av_train_flag = rnd.Next(0, 10);  // creates a number between 0 and 9
                        if (av_train_flag < 2)
                        {
                            trainfile_writer.Write("-1");
                            for (int k = 1; k < 41; ++k)
                            {
                                trainfile_writer.Write("\t" + k + ":" + av_feature_list[m][k]);
                            }
                            trainfile_writer.Write("\n");
                        }
                        else
                        {
                            testfile_av_writer.Write("-1");
                            for (int k = 1; k < 41; ++k)
                            {
                                testfile_av_writer.Write("\t" + k + ":" + av_feature_list[m][k]);
                            }
                            testfile_av_writer.Write("\n");
                        }
                    }
                    
                                       
                    trainfile_writer.Close();
                    testfile_writer.Close();
                    earliestfile_writer.Close();
                    testfile_av_writer.Close();
                }
            }
        }

        static void Main(string[] args) //Main_12_11To1_lineage
        {
            string allfeature_file = "../genomic_features/all_features_seqname_lineage_year.tsv";
            List<String[]> ls = ReadTSV(allfeature_file);
            List<String[]> pe_ls = new List<String[]>();
            List<String[]> sp_ls = new List<String[]>();
            List<String[]> av_ls = new List<String[]>();
            HashSet<string> pe_lineages = new HashSet<string>(new string[] { "Hu1", "Hu2", "Hu3", "Hu4", "Eq", "Ca1", "Ca2", "Sw1", "Sw2", "Sw3", "Sw4", "Sw5" });
            for (int i = 0; i < ls.Count; ++i)
            {
                if (pe_lineages.Contains(ls[i][42])) pe_ls.Add(ls[i]);
                else if (ls[i][42] == "Av") av_ls.Add(ls[i]);
                else sp_ls.Add(ls[i]);
            }
            //pe_lineage2features
            Dictionary<string, List<String[]>> lineage_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < pe_ls.Count; ++i)
            {
                string lineage = pe_ls[i][42];
                if (!lineage_features.ContainsKey(lineage))
                {
                    lineage_features.Add(lineage, new List<String[]>());
                }
                lineage_features[lineage].Add(pe_ls[i]);
            }
            Console.WriteLine(lineage_features.Count);

            //list of lineage features
            List<List<String[]>> lineage_feature_list = new List<List<String[]>>();
            foreach (var item in lineage_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                lineage_feature_list.Add(sort_list);
            }
            Console.WriteLine(lineage_feature_list.Count);

            //sporadic
            //sp2features
            Dictionary<string, List<String[]>> sp_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < sp_ls.Count; ++i)
            {
                string sp = sp_ls[i][42];
                if (!sp_features.ContainsKey(sp))
                {
                    sp_features.Add(sp, new List<String[]>());
                }
                sp_features[sp].Add(sp_ls[i]);
            }
            Console.WriteLine(sp_features.Count);

            //list of sp features
            List<String[]> sp_feature_list = new List<String[]>();
            foreach (var item in sp_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                sp_feature_list.AddRange(sort_list);
            }
            Console.WriteLine(sp_feature_list.Count);

            //avian
            //av2features
            Dictionary<string, List<String[]>> av_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < av_ls.Count; ++i)
            {
                string av = av_ls[i][41].Split('|')[2]; //AvHxNx
                if (!av_features.ContainsKey(av))
                {
                    av_features.Add(av, new List<String[]>());
                }
                av_features[av].Add(av_ls[i]);
            }
            Console.WriteLine(av_features.Count);

            //list of av features
            List<String[]> av_feature_list = new List<String[]>();
            foreach (var item in av_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                av_feature_list.AddRange(sort_list);
            }
            Console.WriteLine(av_feature_list.Count);

            int datasetId = 0;
            Random rnd = new Random();
            for (int lineage_index_1 = 0; lineage_index_1 < lineage_feature_list.Count; lineage_index_1++)
            {
                    string trainfile = "12_11To1_lineage/dataset_12_11To1Lineage/train_dataset_" + datasetId;
                    string testfile = "12_11To1_lineage/dataset_12_11To1Lineage/test_dataset_" + datasetId;
                    string testfile_av = "12_11To1_lineage/dataset_12_11To1Lineage/av_dataset_" + datasetId;
                    StreamWriter trainfile_writer = new StreamWriter(trainfile, false, Encoding.ASCII);
                    StreamWriter testfile_writer = new StreamWriter(testfile, false, Encoding.ASCII);
                    StreamWriter testfile_av_writer = new StreamWriter(testfile_av, false, Encoding.ASCII);
                    datasetId++;
                    //positive
                    for (int j = 0; j < 12; ++j)
                    {
                        if (j != lineage_index_1)
                        {
                            for (int m = 0; m < lineage_feature_list[j].Count; ++m)
                            {
                                trainfile_writer.Write("1");
                                for (int k = 1; k < 41; ++k)
                                {
                                    trainfile_writer.Write("\t" + k + ":" + lineage_feature_list[j][m][k]);
                                }
                                trainfile_writer.Write("\n");
                            }
                        }
                    }
                    for (int m = 0; m < lineage_feature_list[lineage_index_1].Count; ++m)
                    {
                        testfile_writer.Write("1");
                        for (int k = 1; k < 41; ++k)
                        {
                            testfile_writer.Write("\t" + k + ":" + lineage_feature_list[lineage_index_1][m][k]);
                        }
                        testfile_writer.Write("\n");
                    }

                    //negative samples: 80% for training, 20% for testing
                    //sporadic
                    for (int m = 0; m < sp_feature_list.Count; ++m)
                    {
                        int sp_train_flag = rnd.Next(0, 10);  // creates a number between 0 and 9
                        if (sp_train_flag < 8)
                        {
                            trainfile_writer.Write("-1");
                            for (int k = 1; k < 41; ++k)
                            {
                                trainfile_writer.Write("\t" + k + ":" + sp_feature_list[m][k]);
                            }
                            trainfile_writer.Write("\n");
                        }
                        else
                        {
                            testfile_writer.Write("-1");
                            for (int k = 1; k < 41; ++k)
                            {
                                testfile_writer.Write("\t" + k + ":" + sp_feature_list[m][k]);
                            }
                            testfile_writer.Write("\n");
                        }
                    }

                    //avian
                    for (int m = 0; m < av_feature_list.Count; ++m)
                    {
                        int av_train_flag = rnd.Next(0, 10);  // creates a number between 0 and 9
                        if (av_train_flag < 2)
                        {
                            trainfile_writer.Write("-1");
                            for (int k = 1; k < 41; ++k)
                            {
                                trainfile_writer.Write("\t" + k + ":" + av_feature_list[m][k]);
                            }
                            trainfile_writer.Write("\n");
                        }
                        else
                        {
                            testfile_av_writer.Write("-1");
                            for (int k = 1; k < 41; ++k)
                            {
                                testfile_av_writer.Write("\t" + k + ":" + av_feature_list[m][k]);
                            }
                            testfile_av_writer.Write("\n");
                        }
                    }

                    trainfile_writer.Close();
                    testfile_writer.Close();
                    testfile_av_writer.Close();
            }
        }

        static void Main_noHu3Hu4_lineage(string[] args) //Main_noHu3Hu4_lineage
        {
            string allfeature_file = "../genomic_features/all_features_seqname_lineage_year.tsv";
            List<String[]> ls = ReadTSV(allfeature_file);
            List<String[]> pe_ls = new List<String[]>();
            List<String[]> sp_ls = new List<String[]>();
            List<String[]> av_ls = new List<String[]>();
            HashSet<string> pe_lineages = new HashSet<string>(new string[] { "Hu1", "Hu2", "Hu3", "Hu4", "Eq", "Ca1", "Ca2", "Sw1", "Sw2", "Sw3", "Sw4", "Sw5" });
            for (int i = 0; i < ls.Count; ++i)
            {
                if (pe_lineages.Contains(ls[i][42])) pe_ls.Add(ls[i]);
                else if (ls[i][42] == "Av") av_ls.Add(ls[i]);
                else sp_ls.Add(ls[i]);
            }
            //pe_lineage2features
            Dictionary<string, List<String[]>> lineage_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < pe_ls.Count; ++i)
            {
                string lineage = pe_ls[i][42];
                if (!lineage_features.ContainsKey(lineage))
                {
                    lineage_features.Add(lineage, new List<String[]>());
                }
                lineage_features[lineage].Add(pe_ls[i]);
            }
            Console.WriteLine(lineage_features.Count);

            //list of lineage features
            List<List<String[]>> lineage_feature_list = new List<List<String[]>>();
            foreach (var item in lineage_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                lineage_feature_list.Add(sort_list);
            }
            Console.WriteLine(lineage_feature_list.Count);

            //sporadic
            //sp2features
            Dictionary<string, List<String[]>> sp_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < sp_ls.Count; ++i)
            {
                string sp = sp_ls[i][42];
                if (!sp_features.ContainsKey(sp))
                {
                    sp_features.Add(sp, new List<String[]>());
                }
                sp_features[sp].Add(sp_ls[i]);
            }
            Console.WriteLine(sp_features.Count);

            //list of sp features
            List<String[]> sp_feature_list = new List<String[]>();
            foreach (var item in sp_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                sp_feature_list.AddRange(sort_list);
            }
            Console.WriteLine(sp_feature_list.Count);

            //avian
            //av2features
            Dictionary<string, List<String[]>> av_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < av_ls.Count; ++i)
            {
                string av = av_ls[i][41].Split('|')[2]; //AvHxNx
                if (!av_features.ContainsKey(av))
                {
                    av_features.Add(av, new List<String[]>());
                }
                av_features[av].Add(av_ls[i]);
            }
            Console.WriteLine(av_features.Count);

            //list of av features
            List<String[]> av_feature_list = new List<String[]>();
            foreach (var item in av_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                av_feature_list.AddRange(sort_list);
            }
            Console.WriteLine(av_feature_list.Count);

            string earliest_lineage_file = "../genomic_features/earliest_features.tsv";
            List<String[]> earliest_feature_list = ReadTSV(earliest_lineage_file);

            Random rnd = new Random();

            for (int i = 0; i < 1; ++i)
            {
                int lineage_index_1 = 2;  // Hu3
                int lineage_index_2 = 3;  // Hu4
                string trainfile = "noHu3Hu4_lineage/dataset/train_dataset_" + i;
                string testfile = "noHu3Hu4_lineage/dataset/test_dataset_" + i;
                string earliestfile = "noHu3Hu4_lineage/dataset/earliest_dataset_" + i;
                string testfile_av = "noHu3Hu4_lineage/dataset/av_dataset_" + i;
                StreamWriter trainfile_writer = new StreamWriter(trainfile, false, Encoding.ASCII);
                StreamWriter testfile_writer = new StreamWriter(testfile, false, Encoding.ASCII);
                StreamWriter earliestfile_writer = new StreamWriter(earliestfile, false, Encoding.ASCII);
                StreamWriter testfile_av_writer = new StreamWriter(testfile_av, false, Encoding.ASCII);
                //positive
                for (int j = 0; j < 12; ++j)
                {
                    if (j != lineage_index_1 && j != lineage_index_2)
                    {
                        for (int m = 0; m < lineage_feature_list[j].Count; ++m)
                        {
                            trainfile_writer.Write("1");
                            for (int k = 1; k < 41; ++k)
                            {
                                trainfile_writer.Write("\t" + k + ":" + lineage_feature_list[j][m][k]);
                            }
                            trainfile_writer.Write("\n");
                        }
                    }
                }
                for (int m = 0; m < lineage_feature_list[lineage_index_1].Count; ++m)
                {
                    testfile_writer.Write("1");
                    for (int k = 1; k < 41; ++k)
                    {
                        testfile_writer.Write("\t" + k + ":" + lineage_feature_list[lineage_index_1][m][k]);
                    }
                    testfile_writer.Write("\n");
                }
                for (int m = 0; m < lineage_feature_list[lineage_index_2].Count; ++m)
                {
                    testfile_writer.Write("1");
                    for (int k = 1; k < 41; ++k)
                    {
                        testfile_writer.Write("\t" + k + ":" + lineage_feature_list[lineage_index_2][m][k]);
                    }
                    testfile_writer.Write("\n");
                }

                //earliest sequences for test
                earliestfile_writer.Write("1");
                for (int k = 1; k < 41; ++k)
                {
                    earliestfile_writer.Write("\t" + earliest_feature_list[lineage_index_1][k]);
                }
                earliestfile_writer.Write("\n");
                earliestfile_writer.Write("1");
                for (int k = 1; k < 41; ++k)
                {
                    earliestfile_writer.Write("\t" + earliest_feature_list[lineage_index_2][k]);
                }
                earliestfile_writer.Write("\n");

                //negative samples: 80% for training, 20% for testing
                //sporadic
                for (int m = 0; m < sp_feature_list.Count; ++m)
                {
                    int sp_train_flag = rnd.Next(0, 10);  // creates a number between 0 and 9
                    if (sp_train_flag < 8)
                    {
                        trainfile_writer.Write("-1");
                        for (int k = 1; k < 41; ++k)
                        {
                            trainfile_writer.Write("\t" + k + ":" + sp_feature_list[m][k]);
                        }
                        trainfile_writer.Write("\n");
                    }
                    else
                    {
                        testfile_writer.Write("-1");
                        for (int k = 1; k < 41; ++k)
                        {
                            testfile_writer.Write("\t" + k + ":" + sp_feature_list[m][k]);
                        }
                        testfile_writer.Write("\n");
                    }
                }
        
                for (int m = 0; m < av_feature_list.Count; ++m)
                {
                    testfile_av_writer.Write("-1");
                    for (int k = 1; k < 41; ++k)
                    {
                        testfile_av_writer.Write("\t" + k + ":" + av_feature_list[m][k]);
                    }
                    testfile_av_writer.Write("\n");
                }

                trainfile_writer.Close();
                testfile_writer.Close();
                earliestfile_writer.Close();
                testfile_av_writer.Close();
            }
        }
        

        static void Main_training_dataset_80To20(string[] args) //_training_dataset_80To20, main model
        {
             string allfeature_file = "../genomic_features/all_features_seqname_lineage_year.tsv";
            List<String[]> ls = ReadTSV(allfeature_file);
            List<String[]> pe_ls = new List<String[]>();
            List<String[]> sp_ls = new List<String[]>();
            List<String[]> av_ls = new List<String[]>();
            HashSet<string> pe_lineages = new HashSet<string>(new string[] { "Hu1", "Hu2", "Hu3", "Hu4", "Eq", "Ca1", "Ca2", "Sw1", "Sw2", "Sw3", "Sw4", "Sw5" });
            for (int i = 0; i < ls.Count; ++i)
            {
                if (pe_lineages.Contains(ls[i][42])) pe_ls.Add(ls[i]);
                else if (ls[i][42] == "Av") av_ls.Add(ls[i]);
                else sp_ls.Add(ls[i]);
            }
            //pe_lineage2features
            Dictionary<string, List<String[]>> lineage_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < pe_ls.Count; ++i)
            {
                string lineage = pe_ls[i][42];
                if (!lineage_features.ContainsKey(lineage))
                {
                    lineage_features.Add(lineage, new List<String[]>());
                }
                lineage_features[lineage].Add(pe_ls[i]);
            }
            Console.WriteLine(lineage_features.Count);

            //list of lineage features
            List<List<String[]>> lineage_feature_list = new List<List<String[]>>();
            foreach (var item in lineage_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                lineage_feature_list.Add(sort_list);
            }
            Console.WriteLine(lineage_feature_list.Count);

            //sporadic
            //sp2features
            Dictionary<string, List<String[]>> sp_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < sp_ls.Count; ++i)
            {
                string sp = sp_ls[i][42];
                if (!sp_features.ContainsKey(sp))
                {
                    sp_features.Add(sp, new List<String[]>());
                }
                sp_features[sp].Add(sp_ls[i]);
            }
            Console.WriteLine(sp_features.Count);

            //list of sp features
            List<String[]> sp_feature_list = new List<String[]>();
            foreach (var item in sp_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                sp_feature_list.AddRange(sort_list);
            }
            Console.WriteLine(sp_feature_list.Count);

            //avian
            //av2features
            Dictionary<string, List<String[]>> av_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < av_ls.Count; ++i)
            {
                string av = av_ls[i][41].Split('|')[2]; //AvHxNx
                if (!av_features.ContainsKey(av))
                {
                    av_features.Add(av, new List<String[]>());
                }
                av_features[av].Add(av_ls[i]);
            }
            Console.WriteLine(av_features.Count);

            //list of sp features
            List<String[]> av_feature_list = new List<String[]>();
            foreach (var item in av_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                av_feature_list.AddRange(sort_list);
            }
            Console.WriteLine(av_feature_list.Count);

            string earliest_lineage_file = "../genomic_features/earliest_features.tsv";
            List<String[]> earliest_feature_list = ReadTSV(earliest_lineage_file);
            HashSet<string> earliest_virus = new HashSet<string>();
            for(int i = 0; i < earliest_feature_list.Count; ++i)
            {
                earliest_virus.Add(earliest_feature_list[i][41]);
            }

            Random rnd = new Random();

            string trainfile = "/SVM_model/cds/training_features_noEarliest.tsv";
            string testfile = "/SVM_model/cds//testing_features_Earliest.tsv";
            string testfile_av = "/SVM_model/cds/avian_features.tsv";
            StreamWriter trainfile_writer = new StreamWriter(trainfile, false, Encoding.ASCII);
            StreamWriter testfile_writer = new StreamWriter(testfile, false, Encoding.ASCII);
            StreamWriter testfile_av_writer = new StreamWriter(testfile_av, false, Encoding.ASCII);

            //80 % for training, 20 % for testing, excluding the earliest viruses
            //positive
            for (int j = 0; j < 12; ++j)
            {
                for (int m = 0; m < lineage_feature_list[j].Count; ++m)
                {
                    int pe_train_flag = rnd.Next(0, 10);  // creates a number between 0 and 9
                    if (earliest_virus.Contains(lineage_feature_list[j][m][41])) {
                        pe_train_flag = 10;
                        Console.WriteLine(lineage_feature_list[j][m][41]);
                    }
                    if (pe_train_flag < 8)
                    {
                        trainfile_writer.Write("1");
                        for (int k = 1; k < 41; ++k)
                        {
                            trainfile_writer.Write("\t" + k + ":" + lineage_feature_list[j][m][k]);
                        }
                        trainfile_writer.Write("\n");
                    }
                    else
                    {
                        testfile_writer.Write("1");
                        for (int k = 1; k < 41; ++k)
                        {
                            testfile_writer.Write("\t" + k + ":" + lineage_feature_list[j][m][k]);
                        }
                        testfile_writer.Write("\n");
                    }
                } 
            }
                
            //negative samples: 
            //sporadic
            for (int m = 0; m < sp_feature_list.Count; ++m)
            {
                int sp_train_flag = rnd.Next(0, 10);  // creates a number between 0 and 9
                if (sp_train_flag < 8)
                {
                    trainfile_writer.Write("-1");
                    for (int k = 1; k < 41; ++k)
                    {
                        trainfile_writer.Write("\t" + k + ":" + sp_feature_list[m][k]);
                    }
                    trainfile_writer.Write("\n");
                }
                else
                {
                    testfile_writer.Write("-1");
                    for (int k = 1; k < 41; ++k)
                    {
                        testfile_writer.Write("\t" + k + ":" + sp_feature_list[m][k]);
                    }
                    testfile_writer.Write("\n");
                }
            }
         
            for (int m = 0; m < av_feature_list.Count; ++m)
            {
                testfile_av_writer.Write("-1");
                for (int k = 1; k < 41; ++k)
                {
                    testfile_av_writer.Write("\t" + k + ":" + av_feature_list[m][k]);
                }
                testfile_av_writer.Write("\n");
            }

            trainfile_writer.Close();
            testfile_writer.Close();
            testfile_av_writer.Close();
        }

        static void Main_later2early(string[] args) //_early2later or _later2early
        {
            string allfeature_file = "../genomic_features/all_features_seqname_lineage_year.tsv";
            List<String[]> ls = ReadTSV(allfeature_file);
            List<String[]> pe_ls = new List<String[]>();
            List<String[]> sp_ls = new List<String[]>();
            List<String[]> av_ls = new List<String[]>();
            HashSet<string> pe_lineages = new HashSet<string>(new string[] { "Hu1", "Hu2", "Hu3", "Hu4", "Eq", "Ca1", "Ca2", "Sw1", "Sw2", "Sw3", "Sw4", "Sw5" });
            for (int i = 0; i < ls.Count; ++i)
            {
                if (pe_lineages.Contains(ls[i][42])) pe_ls.Add(ls[i]);
                else if (ls[i][42] == "Av") av_ls.Add(ls[i]);
                else sp_ls.Add(ls[i]);
            }
            //pe_lineage2features
            Dictionary<string, List<String[]>> lineage_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < pe_ls.Count; ++i)
            {
                string lineage = pe_ls[i][42];
                if (!lineage_features.ContainsKey(lineage))
                {
                    lineage_features.Add(lineage, new List<String[]>());
                }
                lineage_features[lineage].Add(pe_ls[i]);
            }
            Console.WriteLine(lineage_features.Count);

            //list of lineage features
            List<List<String[]>> lineage_feature_list = new List<List<String[]>>();
            foreach (var item in lineage_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                lineage_feature_list.Add(sort_list);
            }
            Console.WriteLine(lineage_feature_list.Count);

            //sporadic
            //sp2features
            Dictionary<string, List<String[]>> sp_features = new Dictionary<string, List<String[]>>();
            for (int i = 0; i < sp_ls.Count; ++i)
            {
                string sp = sp_ls[i][42];
                if (!sp_features.ContainsKey(sp))
                {
                    sp_features.Add(sp, new List<String[]>());
                }
                sp_features[sp].Add(sp_ls[i]);
            }
            Console.WriteLine(sp_features.Count);

            //list of sp features
            List<String[]> sp_feature_list = new List<String[]>();
            foreach (var item in sp_features)
            {
                List<String[]> sort_list = item.Value;
                sort_list.Sort((x, y) => Double.Parse(x[43]).CompareTo(Double.Parse(y[43])));
                sp_feature_list.AddRange(sort_list);
            }
            Console.WriteLine(sp_feature_list.Count);

            string trainfile = "later2early/dataset/train_dataset";
            string testfile = "later2early/dataset/test_dataset";
            StreamWriter trainfile_writer = new StreamWriter(trainfile, false, Encoding.ASCII);
            StreamWriter testfile_writer = new StreamWriter(testfile, false, Encoding.ASCII);
            Random rnd = new Random();
            //positive
            for (int j = 0; j < 12; ++j)
            {
                for (int m = 0; m < lineage_feature_list[j].Count; ++m)
                {
                    //current is later2early
                    //if use early2later, just exchange the testfile_writer and trainfile_writer, and set the percentage = 0.8
                    if (m < lineage_feature_list[j].Count * 0.2) 
                    {
                        testfile_writer.Write("1");
                        for (int k = 1; k < 41; ++k)
                        {
                            testfile_writer.Write("\t" + k + ":" + lineage_feature_list[j][m][k]);
                        }
                        testfile_writer.Write("\n");
                    }
                    else
                    {
                        trainfile_writer.Write("1");
                        for (int k = 1; k < 41; ++k)
                        {
                            trainfile_writer.Write("\t" + k + ":" + lineage_feature_list[j][m][k]);
                        }
                        trainfile_writer.Write("\n");
                    }
                }
            }

            //negative samples: 80% for training, 20% for testing
            for (int m = 0; m < sp_feature_list.Count; ++m)
            {
                int sp_train_flag = rnd.Next(0, 10);  // creates a number between 0 and 9
                if (sp_train_flag < 8)
                {
                    trainfile_writer.Write("-1");
                    for (int k = 1; k < 41; ++k)
                    {
                        trainfile_writer.Write("\t" + k + ":" + sp_feature_list[m][k]);
                    }
                    trainfile_writer.Write("\n");
                }
                else
                {
                    testfile_writer.Write("-1");
                    for (int k = 1; k < 41; ++k)
                    {
                        testfile_writer.Write("\t" + k + ":" + sp_feature_list[m][k]);
                    }
                    testfile_writer.Write("\n");
                }
            }

            trainfile_writer.Close();
            testfile_writer.Close();
        }

        static List<String[]> ReadTSV(string filePathName)
        {
            List<String[]> ls = new List<String[]>();
            StreamReader fileReader = new StreamReader(filePathName, Encoding.Default);
            string line = "";
            while ((line = fileReader.ReadLine()) != null)
            {
                line = line.Trim();
                line = line.Replace("\"", "");
                if (line.Length > 0)
                {
                    ls.Add(line.Split('\t'));
                }
            }
            fileReader.Close();
            return ls;
        }

        static List<String[]> ReadCSV(string filePathName)
        {
            List<String[]> ls = new List<String[]>();
            StreamReader fileReader = new StreamReader(filePathName, Encoding.Default);
            string line = "";
            while ((line = fileReader.ReadLine()) != null)
            {
                line = line.Trim();
                line = line.Replace("\"", "");
                if (line.Length > 0)
                {
                    ls.Add(line.Split(','));
                }
            }
            fileReader.Close();
            return ls;
        }
    }
}
