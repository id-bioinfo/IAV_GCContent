using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RNAfold_script
{
    class RNAfold_script
    {
        static void Main_RNAfold_extract(string[] args) //_RNAfold_extract
        {
            string rnafold_outfile = args[0]; ///data/ytye/influenza/RNAFold/rnafold/HA_rnafold.txt
            Dictionary<string, string> seqnanme2seq = new Dictionary<string, string>();
            Dictionary<string, List<double>> seqnanme2info = new Dictionary<string, List<double>>();
            int i = 0;
            List<String[]> ls = ReadSpace(rnafold_outfile);
            while (i < ls.Count)
            {
                if (ls[i][0][0] == '>') //header
                {
                    string seqname = ls[i][0];
                    seqnanme2seq.Add(seqname, ls[i + 1][0]);
                    double mfe = Double.Parse(ls[i + 2][1].Substring(1, ls[i + 2][1].Length - 2));
                    double collecttime = Double.Parse(seqname.Split('|')[3]);
                    List<double> temp = new List<double>();
                    temp.Add(collecttime);
                    temp.Add(mfe);
                    seqnanme2info.Add(seqname, temp);
                }
                i++;
            }
            Console.WriteLine(seqnanme2info.Count);
            Parallel.ForEach(seqnanme2seq, item =>
            {
                string seqname = item.Key;
                int[] res = getGCCountAndLength_exact(item.Value);
                double gccontent = res[0] / (double)res[1] * 100;
                seqnanme2info[seqname].Add(gccontent);
                seqnanme2info[seqname].Add(res[1]); //length ATCG
                seqnanme2info[seqname].Add(item.Value.Length);

                res = getGC3CountAndLength_exact(item.Value);
                double gc3 = res[0] / (double)res[1] * 100;
                seqnanme2info[seqname].Add(gc3);

                res = getGC12CountAndLength_exact(item.Value);
                double gc12 = res[0] / (double)res[1] * 100;
                seqnanme2info[seqname].Add(gc12);

            });
            string rnafold_statics_file = rnafold_outfile.Replace(".txt", "") + "_statics.csv";
            StreamWriter fileWriter = new StreamWriter(rnafold_statics_file, false, Encoding.ASCII);
            foreach (var item in seqnanme2info)
            {
                fileWriter.Write(item.Key.Substring(1));
                for (int j = 0; j < item.Value.Count; ++j)
                    fileWriter.Write("," + item.Value[j]);
                fileWriter.Write("\n");
            }
            fileWriter.Close();
        }

        static void Main_linage(string[] args)
        {
            string rnafold_statics_file = args[0];
            List<String[]> ls = ReadCSV(rnafold_statics_file);
            Dictionary<string, String[]> seqnanme2info = new Dictionary<string, String[]>();
            for (int i = 0; i < ls.Count; i++)
            {
                seqnanme2info.Add(ls[i][0], ls[i]);
            }
            Console.WriteLine(ls.Count);
            string lineagefile = args[1];
            List<String[]> gdata = ReadTSV(lineagefile);
            Console.WriteLine(gdata.Count);

            string[] usedlineage = { "Hu1", "Hu2", "Hu3", "Hu4", "Eq", "Ca1", "Ca2", "Sw1", "Sw2", "Sw3", "Sw4", "Sw5" };

            HashSet<string> usedlineageSet = new HashSet<string>(usedlineage);


            Dictionary<string, List<String>> lineageAll = new Dictionary<string, List<String>>();
            for (int i = 0; i < gdata.Count; ++i)
            {
                string lineage = gdata[i][0];
                if (usedlineageSet.Contains(lineage))
                {
                    lineageAll.Add(lineage, new List<String>());
                }
               
                for (int j = 1; j < gdata[i].Length; ++j)
                {
                    lineageAll[lineage].Add(gdata[i][j]);
                }
            }
            string outfilePrefix = rnafold_statics_file.Replace(".txt", "");
            outfilePrefix = outfilePrefix.Replace(".csv", "");
            outfilePrefix = "lineage/" + outfilePrefix;
            foreach (var item in lineageAll)
            {
                string lineage = item.Key;
                string outfile = outfilePrefix + "_" + lineage + ".csv";
                StreamWriter fileWriter = new StreamWriter(outfile, false, Encoding.ASCII);
                for (int j = 0; j < item.Value.Count; ++j)
                {
                    string seqname = item.Value[j];
                    if (!seqnanme2info.ContainsKey(seqname))
                    {
                        Console.WriteLine("error can not find " + seqname + " in " + rnafold_statics_file);
                        continue;
                    }
                    fileWriter.WriteLine(string.Join(",", seqnanme2info[seqname]));
                }
                fileWriter.Close();
            }
        }

        static void Main_summary(string[] args)
        {
            string rnafold_statics_file = args[0];
            List<String[]> ls = ReadCSV(rnafold_statics_file);
            Dictionary<string, String[]> seqnanme2info = new Dictionary<string, String[]>();
            for (int i = 0; i < ls.Count; i++)
            {
                seqnanme2info.Add(ls[i][0], ls[i]);
            }
            Console.WriteLine(ls.Count);
            string lineagefile = args[1];
            List<String[]> gdata = ReadTSV(lineagefile);
            Console.WriteLine(gdata.Count);

            string[] mammal_lineage = { "Hu1", "Hu2", "Hu3", "Hu4", "Eq", "Ca1", "Ca2", "Sw1", "Sw2", "Sw3", "Sw4", "Sw5" };

            HashSet<string> mammallineageSet = new HashSet<string>(mammal_lineage);

            Dictionary<string, List<String>> lineageAll = new Dictionary<string, List<String>>();
            for (int i = 0; i < gdata.Count; ++i)
            {
                string lineage = gdata[i][0];
                if (mammallineageSet.Contains(lineage)) lineage = "Ma";
                else if (lineage != "Av") lineage = "Sp";
                if (!lineageAll.ContainsKey(lineage)) lineageAll.Add(lineage, new List<String>());
                for (int j = 1; j < gdata[i].Length; ++j)
                {
                    lineageAll[lineage].Add(gdata[i][j]);
                }
            }
            string outfilePrefix = rnafold_statics_file.Replace(".txt", "");
            outfilePrefix = outfilePrefix.Replace(".csv", "");
            outfilePrefix = "lineage/" + outfilePrefix;
            string outfile = outfilePrefix + "_MFE_GC.csv";
            StreamWriter fileWriter = new StreamWriter(outfile, false, Encoding.ASCII);
            fileWriter.WriteLine("lineage,seqname,collectiontime,mfe,gc,atcgSeqL,seqL");
            foreach (var item in lineageAll)
            {
                string lineage = item.Key;
                for (int j = 0; j < item.Value.Count; ++j)
                {
                    string seqname = item.Value[j];
                    if (!seqnanme2info.ContainsKey(seqname))
                    {
                        Console.WriteLine("error can not find " + seqname + " in " + rnafold_statics_file);
                        continue;
                    }
                    fileWriter.WriteLine(lineage + "," + string.Join(",", seqnanme2info[seqname]));
                }
            }
            fileWriter.Close();
        }

        //add GC3 and GC12
        static void Main_add_GC3_GC12(string[] args) //add_GC3_GC12
        {
            string gene = args[0]; //HA
            //read msa
            string msafile = "/data/ytye/influenza/RNAFold/pal2nal/" + gene + "_nt_all.fasta";
            Dictionary<string, string> seqnanme2seq = ReadMSA(msafile);

            //read lineage, mfe and others
            string[] lineages = { "Av", "Ca1", "Ca2", "Eq", "Hu1", "Hu2", "Hu3", "Hu4", "Sp", "Sw1", "Sw2", "Sw3", "Sw4", "Sw5" };
            string prefixPath = "/data/ytye/influenza/RNAFold/lineage/";
            for (int i = 0; i < lineages.Length; ++i)
            {
                string mfefile = prefixPath + gene + "_rnafold_statics_" + lineages[i] + ".csv";
                List<String[]> ls = ReadCSV(mfefile);
                string rnafold_statics_file = prefixPath + gene + "_rnafold_addGC3_12_statics_" + lineages[i] + ".csv";
                StreamWriter fileWriter = new StreamWriter(rnafold_statics_file, false, Encoding.ASCII);

                for (int j = 0; j < ls.Count; ++j)
                {
                    string seqname = ">" + ls[j][0];
                    if (!seqnanme2seq.ContainsKey(seqname))
                        Console.WriteLine("can not find " + seqname + " in " + gene);
                    string seq = seqnanme2seq[seqname];
                    int[] res = getGC3CountAndLength_exact(seq);
                    double gc3 = res[0] / (double)res[1] * 100;
                    res = getGC12CountAndLength_exact(seq);
                    double gc12 = res[0] / (double)res[1] * 100;
                    //res = getGCCountAndLength_exact(seq);
                    //double gc = res[0] / (double)res[1] * 100;
                    fileWriter.WriteLine(string.Join(",", ls[j]) + "," + gc3 + "," + gc12);
                }
                fileWriter.Close();
            }
        }

        static void Main_lineageSeq(string[] args) //get lineage sequences
        {
            string[] genes = { "PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1" };
            Dictionary<string, Dictionary<string, string>> allseqnanme2seq = new Dictionary<string, Dictionary<string, string>>();
            for (int k = 0; k < 8; ++k)
            {
                string gene = genes[k];
                string msafile = "/data/ytye/influenza/RNAFold/pal2nal/" + gene + "_nt_all.fasta";
                Dictionary<string, string> seqnanme2seq = ReadMSA(msafile);
                allseqnanme2seq.Add(gene, seqnanme2seq);
            }

            //read lineage, mfe and others
            string[] lineages = { "Ca1", "Ca2", "Eq", "Hu1", "Hu2", "Hu3", "Hu4", "Sw1", "Sw2", "Sw3", "Sw4", "Sw5"};
            string prefixPath = "/data/ytye/influenza/RNAFold/lineage/";

            for (int i = 0; i < lineages.Length; ++i)
            {
                Dictionary<string, string> sequences = new Dictionary<string, string>();
                Dictionary<string, int> sequence2gene = new Dictionary<string, int>();
                string mfefile = prefixPath + "HA_rnafold_statics_" + lineages[i] + ".csv";
                List<String[]> ls = ReadCSV(mfefile);
                for (int k = 0; k < 8; ++k)
                {
                    string gene = genes[k];
                    Dictionary<string, string> seqnanme2seq = allseqnanme2seq[gene];
                    for (int j = 0; j < ls.Count; ++j)
                    {
                        string seqname = ">" + ls[j][0];
                        string seq = seqnanme2seq[seqname];
                        if (!sequences.ContainsKey(seqname))
                        {
                            sequences.Add(seqname, "");
                            sequence2gene.Add(seqname, 0);
                        }
                        sequences[seqname] += seq;
                        sequence2gene[seqname]++;
                    }
                }
                string outfile = "/data/ytye/influenza/codon_bias/lineage_seq/" + lineages[i] + "_genome.fas";
                StreamWriter fileWriter = new StreamWriter(outfile, false, Encoding.ASCII);
                foreach (var item in sequences)
                {
                    if (sequence2gene[item.Key] == 8)
                    {
                        fileWriter.WriteLine(item.Key);
                        fileWriter.WriteLine(item.Value);
                    }
                    else
                    {
                        Console.WriteLine("can not find " + item.Key);
                    }
                }
                fileWriter.Close();
            }
        }


        static void Main(string[] args) //_split_time
        { 
            //read lineage, mfe and others
            string[] lineages = { "Ca1", "Ca2", "Eq", "Hu1", "Hu2", "Hu3", "Hu4", "Sw1", "Sw2", "Sw3", "Sw4", "Sw5" };
            string prefixPath = @"D:\work\projects\Influenza evolution\codon_bias\";
            for (int i = 0; i < lineages.Length; ++i)
            {
                string caifile = prefixPath +  "cai_" + lineages[i] + ".csv";
                List<String[]> ls = ReadCSV(caifile);
                string outfile = prefixPath + lineages[i] + "_cai_time" + ".csv";

                StreamWriter fileWriter = new StreamWriter(outfile, false, Encoding.ASCII);

                for (int j = 0; j < ls.Count; ++j)
                {
                    string seqname =ls[j][0];
                    string[] arr = seqname.Split('|');
                    fileWriter.WriteLine(string.Join(",", ls[j]) + "," + arr[3]);
                }
                fileWriter.Close();
            }
        }


        static int[] getGCCountAndLength_exact(string seq)
        {
            int count = 0;
            int totalL = 0;
            int IUPAC = 0;
            for (int i = 0; i < seq.Length; ++i)
            {
                if (seq[i] == 'G' || seq[i] == 'C') count++;
                if (seq[i] == 'G' || seq[i] == 'C' || seq[i] == 'A' || seq[i] == 'T' || seq[i] == 'U') totalL++; //modificy on 20230329
                else IUPAC++;
            }
            int[] res = new int[3];
            res[0] = count;
            res[1] = totalL;
            res[2] = IUPAC;
            return res;
        }

        static int[] getGC3CountAndLength_exact(string seq, bool excludeGap = true)
        {
            int count = 0;
            int totalL = 0;
            if (seq.Length % 3 != 0)
            {
                Console.WriteLine("seqlength error = " + seq.Length + "," + seq);
            }
            for (int i = 2; i < seq.Length; i += 3)
            {
                if (seq[i] == 'G' || seq[i] == 'C') count++;
                if (seq[i] == 'G' || seq[i] == 'C' || seq[i] == 'A' || seq[i] == 'T' || seq[i] == 'U') totalL++;
            }
            int[] res = new int[2];
            res[0] = count;
            res[1] = totalL;
            return res;
        }

        static int[] getGC12CountAndLength_exact(string seq)
        {
            int count = 0;
            int totalL = 0;
            if (seq.Length % 3 != 0)
            {
                Console.WriteLine("seqlength error = " + seq.Length + "," + seq);
            }
            for (int i = 0; i < seq.Length; i++)
            {
                if (i % 3 == 2) continue;
                if (seq[i] == 'G' || seq[i] == 'C') count++;
                if (seq[i] == 'G' || seq[i] == 'C' || seq[i] == 'A' || seq[i] == 'T' || seq[i] == 'U') totalL++;
            }
            int[] res = new int[2];
            res[0] = count;
            res[1] = totalL;
            return res;
        }

        static List<String[]> ReadSpace(string filePathName)
        {
            List<String[]> ls = new List<String[]>();
            StreamReader fileReader = new StreamReader(filePathName, Encoding.Default);
            string line = "";
            while ((line = fileReader.ReadLine()) != null)
            {
                line = line.Trim();
                if (line.Length > 0)
                {
                    ls.Add(line.Split(' '));
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

        static Dictionary<string, string> ReadMSA(string msafile, int minLength = 1, int maxLength = 10000000)
        {
            Dictionary<string, string> msaMap = new Dictionary<string, string>();
            StreamReader fileReader = new StreamReader(msafile, Encoding.Default);
            string line, seqName = "", seqLine = "";
            while ((line = fileReader.ReadLine()) != null)
            {
                line = line.Trim();
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
