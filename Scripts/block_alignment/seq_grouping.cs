//this is a C# program
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace seq_grouping
{
    class seq_grouping
    {
        static void Main(string[] args)
        {
            List<string> proteinId = new List<string> { "HA", "NA", "NP", "PA", "PB1", "PB2", "M1", "NS1" };
            string prefixPath = "filter/";
            string prefixoutPath = "filter/partition/";
            string seqName, subType;
            double collectDate = 0;
            for(int i=0; i<proteinId.Count; ++i)
            {
                string msafile = prefixPath + proteinId[i] + "_aa.fasta";
                Dictionary<string, string> msa = ReadMSA(msafile);
                Console.WriteLine("ntMsa = " + msa.Count);
                string FileUrl = prefixoutPath + "/" + proteinId[i] + "/";
                Directory.CreateDirectory(FileUrl);

                //H3: NorthAmerica_1 (0), NorthAmerica_2 (1), NorthAmerica_3 (2), Asia (3), Europe (4), Other (5)
                //H1: NorthAmerica_1 (6), NorthAmerica_2 (7), Asia (8), Europe (9), Other (10)
                //H5 (11)
                //H9 (12)
                //other (13)
                //20220403
                //0 H3N2: sort by time and split to three 46364 (0,1,2)
                //1 H1N1 and time >=2009:  sort by time and split to two 34118 (3,4)
                //2 H1N2: 2606, 5
                //3 H3N8: 1608, 6
                //4 H5N1: 3937, 7
                //5 H5N6: 1629, 8
                //6 H5N8: 2295, 9
                //7 H4N6: 1485, 10
                //8 H9N2: 4070, 11
                //9 H7N9: 1969, 12
                //10 other: 10383, 13

                List<List<string>> subMSAKey = new List<List<string>>();
                for (int j = 0; j < 11; ++j) subMSAKey.Add(new List<string>());
                Dictionary<string, double> h3n2Keys = new Dictionary<string, double>();
                Dictionary<string, double> h1n1Keys = new Dictionary<string, double>();
                Console.WriteLine(subMSAKey.Count);
                foreach (var item in msa)
                {
                    seqName = item.Key;
                    string[] arr = seqName.Split('|');
                    subType = arr[2];
                    collectDate = Double.Parse(arr[3]);
                    switch (subType)
                    {
                        case "H3N2":
                            //subMSAKey[0].Add(seqName);
                            h3n2Keys.Add(seqName, collectDate);
                            break;
                        case "H1N1":
                            //subMSAKey[1].Add(seqName);
                            h1n1Keys.Add(seqName, collectDate);
                            break;
                        case "H1N2":
                            subMSAKey[2].Add(seqName);
                            break;
                        case "H3N8":
                            subMSAKey[3].Add(seqName);
                            break;
                        case "H5N1":
                            subMSAKey[4].Add(seqName);
                            break;
                        case "H5N6":
                            subMSAKey[5].Add(seqName);
                            break;
                        case "H5N8":
                            subMSAKey[6].Add(seqName);
                            break;
                        case "H4N6":
                            subMSAKey[7].Add(seqName);
                            break;
                        case "H9N2":
                            subMSAKey[8].Add(seqName);
                            break;
                        case "H7N9":
                            subMSAKey[9].Add(seqName);
                            break;
                        default: 
                            subMSAKey[10].Add(seqName);
                            break;
                    }
                }
                int outIdx = 0;
                string outputfile_aaPrefix = FileUrl + proteinId[i] + "_aa_";
                //H3N2
                h3n2Keys = h3n2Keys.OrderBy(u => u.Value).ToDictionary(z => z.Key, y => y.Value);
                List<string> keyList = h3n2Keys.Keys.ToList();
                int onesize = keyList.Count / 3;
                StreamWriter fileWriter_aa = new StreamWriter(outputfile_aaPrefix + outIdx + ".fasta", false, Encoding.ASCII);
                for (int j=0; j < onesize; ++j)
                {
                    fileWriter_aa.WriteLine(keyList[j]);
                    fileWriter_aa.WriteLine(msa[keyList[j]]);
                }
                fileWriter_aa.Close();
                outIdx++;
                fileWriter_aa = new StreamWriter(outputfile_aaPrefix + outIdx + ".fasta", false, Encoding.ASCII);
                for (int j = onesize; j < 2*onesize; ++j)
                {
                    fileWriter_aa.WriteLine(keyList[j]);
                    fileWriter_aa.WriteLine(msa[keyList[j]]);
                }
                fileWriter_aa.Close();
                outIdx++;
                fileWriter_aa = new StreamWriter(outputfile_aaPrefix + outIdx + ".fasta", false, Encoding.ASCII);
                for (int j = 2*onesize; j < keyList.Count; ++j)
                {
                    fileWriter_aa.WriteLine(keyList[j]);
                    fileWriter_aa.WriteLine(msa[keyList[j]]);
                }
                fileWriter_aa.Close();

                //H1N1
                h1n1Keys = h1n1Keys.OrderBy(u => u.Value).ToDictionary(z => z.Key, y => y.Value);
                keyList = h1n1Keys.Keys.ToList();
                onesize = keyList.Count / 2;
                outIdx++;
                fileWriter_aa = new StreamWriter(outputfile_aaPrefix + outIdx + ".fasta", false, Encoding.ASCII);
                for (int j = 0; j < onesize; ++j)
                {
                    fileWriter_aa.WriteLine(keyList[j]);
                    fileWriter_aa.WriteLine(msa[keyList[j]]);
                }
                fileWriter_aa.Close();
                outIdx++;
                fileWriter_aa = new StreamWriter(outputfile_aaPrefix + outIdx + ".fasta", false, Encoding.ASCII);
                for (int j = onesize; j < keyList.Count; ++j)
                {
                    fileWriter_aa.WriteLine(keyList[j]);
                    fileWriter_aa.WriteLine(msa[keyList[j]]);
                }
                fileWriter_aa.Close();

                //2-10
                for (int k = 2; k < 11; ++k)
                {
                    outIdx++;
                    fileWriter_aa = new StreamWriter(outputfile_aaPrefix + outIdx + ".fasta", false, Encoding.ASCII);
                    for (int j = 0; j < subMSAKey[k].Count; ++j)
                    {
                        fileWriter_aa.WriteLine(subMSAKey[k][j]);
                        fileWriter_aa.WriteLine(msa[subMSAKey[k][j]]);
                    }
                    fileWriter_aa.Close();
                }
            }
        }
        
        static Dictionary<string, string> ReadMSA(string msafile, int minLength = 1, int maxLength = 10000000)
        {
            Dictionary<string, string> msaMap = new Dictionary<string, string>();
            StreamReader fileReader = new StreamReader(msafile, Encoding.Default);
            string line, seqName = "", seqLine = "";
            while ((line = fileReader.ReadLine()) != null)
            {
                line = line.Trim();
                if (String.IsNullOrEmpty(line)) continue;
                if (line[0] == '>')
                {
                    if (seqLine.Length > minLength && seqLine.Length < maxLength && !msaMap.ContainsKey(seqName))
                    {
                        msaMap.Add(seqName, seqLine.ToUpper());
                    }
                    seqName = line;
                    seqLine = "";
                }
                else
                {
                    seqLine += line;
                }
            }
            if (seqLine.Length > minLength && seqLine.Length < maxLength && !msaMap.ContainsKey(seqName))
            {
                msaMap.Add(seqName, seqLine.ToUpper());
            }
            fileReader.Close();
            return msaMap;
        }
    }
}
