using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace extractConsensus
{
    class extractConsensus
    {
        static void Main(string[] args)
        {
            List<string> proteinId = new List<string> { "HA", "NA", "NP", "PA", "PB1", "PB2", "M1", "NS1" };
            List<string> consensusName = new List<string> { "H3N2_Sub1", "H3N2_Sub2", "H3N2_Sub3", "H1N1_Sub1", "H1N1_Sub2",
                "H1N2", "H3N8", "H5N1", "H5N6", "H5N8","H4N6","H9N2","H7N9", "Others" };
      
            for (int i = 0; i < proteinId.Count; ++i)
            {
                string prefixPath = "filter/partition/" + proteinId[i] + "/";
                string outFile = prefixPath + proteinId[i] + "_consensus.fasta";
                StreamWriter fileWriter_aa = new StreamWriter(outFile, false, Encoding.ASCII);
                for (int j = 0; j < 14; ++j)
                {
                    string msafile = prefixPath + proteinId[i] + "_aa_all_msa_" + j + ".fasta";
                    Dictionary<string, string> msa = ReadMSA(msafile);
                    string consensus = ExtractConsensusSeq(msa);
                    fileWriter_aa.WriteLine(">" + consensusName[j]);
                    fileWriter_aa.WriteLine(consensus);
                }
                fileWriter_aa.Close();
            }
        }

        static string ExtractConsensusSeq(Dictionary<string, string> msa)
        {
            Dictionary<int, Dictionary<char, int>> pos2CharCount = new Dictionary<int, Dictionary<char, int>>();
            string seqStr;
            foreach (var item in msa)
            {
                seqStr = item.Value;
                for(int i=0; i<seqStr.Length; ++i)
                {
                    if (!pos2CharCount.ContainsKey(i)) pos2CharCount.Add(i, new Dictionary<char, int>());
                    if (!pos2CharCount[i].ContainsKey(seqStr[i])) pos2CharCount[i].Add(seqStr[i], 1);
                    else pos2CharCount[i][seqStr[i]]++;
                }
            }
            char[] seqchars  = msa.First().Value.ToCharArray();
            foreach(var item in pos2CharCount)
            {
                seqchars[item.Key] = getMostFrequencyChar(item.Value); 
            }
            return new String(seqchars);
        }

        static char getMostFrequencyChar(Dictionary<char, int> charCount)
        {
            int maxCount = 0;
            char maxFrequencyChar = '-';
            foreach(var item in charCount)
            {
                if(item.Value > maxCount)
                {
                    maxCount = item.Value;
                    maxFrequencyChar = item.Key;
                }
            }
            return maxFrequencyChar;
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
