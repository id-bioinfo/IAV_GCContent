//This is a C# program.
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace map_aa_nt_meta
{
    class map_aa_nt_meta
    {
        //all potein sequences of 8 poteins are stored in aafile
        //all dna sequences of 8 sgements are stored in nt file
        //metadata should be filtered by laboratory-derived, environmental, and unknown host or subtype sequences before running this program 
        static void Main(string[] args)
        {
            if (args.Length < 4)
            {
                Console.WriteLine("input are aafile ntfile metafile outputdir");
                return;
            }

            List<string> proteinId = new List<string> { "HA", "NA", "NP", "PA", "PB1", "PB2", "M1", "NS1"};
            List<string> geneId = new List<string> { "HA", "NA", "NP", "PA", "PB1", "PB2", "MP", "NS" };
            List<int> proteinLength = new List<int> { 568, 469, 498, 716, 757, 759, 252, 230 };
            double min_alpha = 0.9, max_alpha = 1.1;
            Dictionary<string, double[]> proteinMap = new Dictionary<string, double[]>();
            for (int i = 0; i < proteinId.Count; ++i)
            {
                double[] lengths = new double[2];
                lengths[0] = proteinLength[i] * min_alpha;
                lengths[1] = proteinLength[i] * max_alpha;
                proteinMap.Add(proteinId[i], lengths);
            }
            Dictionary<string, string> geneMap = new Dictionary<string, string>();
            for (int i = 0; i < geneId.Count; ++i)
            {
                geneMap.Add(proteinId[i], geneId[i]);
            }

            //metadata
            string metadatafile = args[2];
            List<String[]> metadataList = ReadCSV(metadatafile);
            Dictionary<string, String[]> metadataMap = new Dictionary<string, String[]>();
            for (int i = 0; i < metadataList.Count; ++i)
            {
                metadataMap.Add(metadataList[i][0], metadataList[i]);
            }
            Console.WriteLine("metadataMap = " + metadataMap.Count);

            //aa
            string aafile = args[0];
            Dictionary<string, string> aaMsa = ReadMSA(aafile);
            Console.WriteLine("aaMsa = " + aaMsa.Count);
            Dictionary<string, Dictionary<string, string>> aaProteinMsa = new Dictionary<string, Dictionary<string, string>>();
            string seqId, aaId, seqName, ntId, seq;
            foreach (var item in aaMsa)
            {
                string[] arr = item.Key.Split('|');
                seqId = arr[0].Substring(1); //remove ">"
                aaId = arr[1];
                if (!metadataMap.ContainsKey(seqId)) continue;
                if (!proteinMap.ContainsKey(aaId)) continue;
                if (item.Value.Length >= proteinMap[aaId][0] && item.Value.Length <= proteinMap[aaId][1])
                {
                    if (!aaProteinMsa.ContainsKey(aaId))
                    {
                        Dictionary<string, string> msas = new Dictionary<string, string>();
                        aaProteinMsa.Add(aaId, msas);
                    }
                    if(!aaProteinMsa[aaId].ContainsKey(seqId))
                        aaProteinMsa[aaId].Add(seqId, item.Value);
                }
            }
         

            //nt
            string ntfile = args[1];
            Dictionary<string, string> ntMsa = ReadMSA(ntfile);
            Console.WriteLine("ntMsa = " + ntMsa.Count);
            
            Dictionary<string, Dictionary<string, string>> ntGeneMsa = new Dictionary<string, Dictionary<string, string>>();
            foreach (var item in ntMsa)
            {
                string[] arr = item.Key.Split('|');
                seqId = arr[0].Substring(1);//remove ">"
                ntId = arr[1];
                if (!ntGeneMsa.ContainsKey(ntId))
                {
                    Dictionary<string, string> msas = new Dictionary<string, string>();
                    ntGeneMsa.Add(ntId, msas);
                }
                if (!ntGeneMsa[ntId].ContainsKey(seqId))
                    ntGeneMsa[ntId].Add(seqId, item.Value);
            }
            //remove identical genome (8 segments)
            Dictionary<string, string> seqId2genome = new Dictionary<string, string>();
            foreach(var item in ntGeneMsa)
            {
                ntId = item.Key;
                foreach(var item1 in item.Value)
                {
                    seqId = item1.Key;
                    seq = item1.Value;
                    if (!seqId2genome.ContainsKey(seqId)) seqId2genome.Add(seqId, "");
                    seqId2genome[seqId] += seq; //combine 8 segments
                }
            }
            List<string> removeSeqIdList = new List<string>();
            HashSet<string> uniqueGenome = new HashSet<string>();
            foreach(var item in seqId2genome)
            {
                if (uniqueGenome.Contains(item.Value))
                    removeSeqIdList.Add(item.Key);
                else
                    uniqueGenome.Add(item.Value);
            }
            foreach (var item in ntGeneMsa)
            {
                ntId = item.Key;
                for (int i = 0; i < removeSeqIdList.Count; ++i)
                {
                    if(ntGeneMsa[ntId].ContainsKey(removeSeqIdList[i]))
                        ntGeneMsa[ntId].Remove(removeSeqIdList[i]);
                }
            }

            //merge
            foreach (var item in aaProteinMsa)
            {
                aaId = item.Key;
                ntId = geneMap[aaId];
                Dictionary<string, string> msas_aa = item.Value;
                Dictionary<string, string> msas_nt = ntGeneMsa[ntId];
                string outputfile_aa = args[3] + aaId + "_aa.fasta";
                string outputfile_nt = args[3] + aaId + "_nt.fasta";
                StreamWriter fileWriter_aa = new StreamWriter(outputfile_aa, false, Encoding.ASCII);
                StreamWriter fileWriter_nt = new StreamWriter(outputfile_nt, false, Encoding.ASCII);
                int k = 0;
                foreach (var item1 in msas_aa)
                {
                    seqId = item1.Key;
                    if (msas_nt.ContainsKey(seqId))
                    {
                        //{ "0:Isolate_Id", "1:Host", "2:Subtype","3:CollectionDate","4:Location","5:Lineage"});
                        seqName = ">" + seqId + "|" + metadataMap[seqId][1] + "|" + metadataMap[seqId][2] + "|" + metadataMap[seqId][3] + "|" + metadataMap[seqId][4].Split('/')[0] + "|" + metadataMap[seqId][5];
                        fileWriter_aa.WriteLine(seqName);
                        fileWriter_aa.WriteLine(item1.Value);
                        fileWriter_nt.WriteLine(seqName);
                        fileWriter_nt.WriteLine(msas_nt[seqId]);
                        k++;
                    }
                }
                Console.WriteLine(aaId + " = " + k);
                fileWriter_aa.Close();
                fileWriter_nt.Close();
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
                if(line[0] == '>')
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
            if(seqLine.Length > minLength && seqLine.Length < maxLength && !msaMap.ContainsKey(seqName))
            {
                msaMap.Add(seqName, seqLine.ToUpper());
            }
            fileReader.Close();
            return msaMap;
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
    }
}
