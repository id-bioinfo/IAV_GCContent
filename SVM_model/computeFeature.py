#install biopython package, https://biopython.org/wiki/Download
from Bio import SeqIO
from collections import defaultdict

def calculate_gc_content(sequence):
    """Calculate GC content of a given DNA sequence."""
    gc_count = 0
    total_count = 0
    for i in range(len(sequence)):
        if sequence[i] == 'G' or sequence[i] == 'C':
            gc_count+=1
        if sequence[i] == 'G' or sequence[i] == 'C' or sequence[i] == 'A' or sequence[i] == 'T':
            total_count+=1

    if total_count == 0:
        return None  # Avoid division by zero
    return gc_count / total_count * 100, gc_count, total_count

def calculate_gc3(sequence):
    if len(sequence)%3 != 0:
        raise ValueError("The length of the sequence is not a multiple of 3.")
    third_positions = sequence[2::3]
    gc_count = 0
    total_count = 0
    #sequence = third_positions
    for i in range(len(sequence)):
        if i%3 == 2:
            if sequence[i] == 'G' or sequence[i] == 'C':
                gc_count+=1
            if sequence[i] == 'G' or sequence[i] == 'C' or sequence[i] == 'A' or sequence[i] == 'T':
                total_count+=1

    if total_count == 0:
        return None  # Avoid division by zero
    return gc_count / total_count * 100, gc_count, total_count

def calculate_gc12(sequence):
    if len(sequence)%3 != 0:
        raise ValueError("The length of the sequence is not a multiple of 3.")
    positions_1_and_2 = sequence[0::3] + sequence[1::3]  # Extract nucleotides at positions 1 and 2

    gc_count = 0
    total_count = 0
    #sequence = positions_1_and_2
    for i in range(len(sequence)):
        if i%3 != 2:
            if sequence[i] == 'G' or sequence[i] == 'C':
                gc_count+=1
            if sequence[i] == 'G' or sequence[i] == 'C' or sequence[i] == 'A' or sequence[i] == 'T':
                total_count+=1

    if total_count == 0:
        return None  # Avoid division by zero
    return gc_count / total_count * 100, gc_count, total_count

def calculate_base_frequency(sequence, base):
    gc_count = 0
    total_count = 0
    for i in range(len(sequence)):
        if sequence[i] == base:
            gc_count+=1
        if sequence[i] == 'G' or sequence[i] == 'C' or sequence[i] == 'A' or sequence[i] == 'T':
            total_count+=1

    if total_count == 0:
        return None  # Avoid division by zero
    return gc_count / total_count * 100, gc_count, total_count

def calculate_dinucleotide_frequency(sequence, first_base, second_base):
    """Calculate dinucleotide frequency of a given DNA sequence."""
    dinucleotide_count = 0
    total_count = 0
    for i in range(len(sequence) - 1):
        if sequence[i] == first_base and sequence[i+1] == second_base:
            dinucleotide_count+=1
        if (sequence[i] == 'G' or sequence[i] == 'C' or sequence[i] == 'A' or sequence[i] == 'T') and \
            (sequence[i+1] == 'G' or sequence[i+1] == 'C' or sequence[i+1] == 'A' or sequence[i+1] == 'T'):
            total_count+=1

    if total_count == 0:
        return 0  # Avoid division by zero
    return dinucleotide_count / total_count * 1000, dinucleotide_count, total_count

# List of 8 FASTA files, one for each segment of the influenza virus
fasta_files = [
    "HA_nt_mink.pal2nal",
    "NA_nt_mink.pal2nal",
    "NP_nt_mink.pal2nal",
    "PA_nt_mink.pal2nal",
    "PB1_nt_mink.pal2nal",
    "PB2_nt_mink.pal2nal",
    "M1_nt_mink.pal2nal",
    "NS1_nt_mink.pal2nal"
]

# Dictionary to store GC content and dinucleotide frequencies for each virus and segment
virus_segment_data = defaultdict(lambda: [None] * len(fasta_files))
virus_genomicGC = defaultdict(lambda: [None] * 2)
virus_genomicGC3 = defaultdict(lambda: [None] * 2)
virus_genomicGC12 = defaultdict(lambda: [None] * 2)
virus_genomicCpG = defaultdict(lambda: [None] * 2)
virus_genomicGpC = defaultdict(lambda: [None] * 2)
virus_genomicGpG = defaultdict(lambda: [None] * 2)
virus_genomicCpC = defaultdict(lambda: [None] * 2)

prefix_path = "demo_minkH5_cds/"
# Process each FASTA file
for i, fasta_file in enumerate(fasta_files, start=0):
    try:
        # Read the FASTA file
        with open(prefix_path + fasta_file, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                # Ensure the sequence is in uppercase
                sequence = str(record.seq).upper()  
                # Calculate GC content for this sequence
                gc_content, gc_count, gc_total_count  = calculate_gc_content(sequence)
                # Calculate CpG frequency for this sequence
                cpg_frequency, cpg_count, cpg_total_count  = calculate_dinucleotide_frequency(sequence, 'C', 'G')
                # Calculate GpC frequency for this sequence
                gpc_frequency, gpc_count, gpc_total_count = calculate_dinucleotide_frequency(sequence, 'G', 'C')
                # Calculate GpG frequency for this sequence
                gpg_frequency, gpg_count, gpg_total_count = calculate_dinucleotide_frequency(sequence, 'G', 'G')
                # Calculate CpC frequency for this sequence
                cpc_frequency, cpc_count, cpc_total_count = calculate_dinucleotide_frequency(sequence, 'C', 'C')

                ##GC3 and GC12
                gc3_content, gc3_count, gc3_total_count = calculate_gc3(sequence)
                gc12_content, gc12_count, gc12_total_count = calculate_gc12(sequence)

                ##O/E ratio
                #c_frequency, c_count, c_total_count = calculate_base_frequency(sequence, 'C')
                #g_frequency, g_count, g_total_count = calculate_base_frequency(sequence, 'C')
                #cpg_oe = cpg_frequency / (c_frequency*g_frequency) * 10 #cpg_frequency is 1000%, c_frequency is 100%, g_frequency is 100%
                #gpc_oe = gpc_frequency / (c_frequency*g_frequency) * 10
                #gpg_oe = gpg_frequency / (g_frequency*g_frequency) * 10
                #cpc_oe = cpc_frequency / (c_frequency*c_frequency) * 10

                # Store the data for this virus and segment (index i)
                virus_segment_data[record.id][i] = (gc_content, cpg_frequency, gpc_frequency, gpg_frequency, cpc_frequency)

                if i == 0:
                    virus_genomicGC[record.id][0] = gc_count
                    virus_genomicGC[record.id][1] = gc_total_count
                    virus_genomicGC3[record.id][0] = gc3_count
                    virus_genomicGC3[record.id][1] = gc3_total_count
                    virus_genomicGC12[record.id][0] = gc12_count
                    virus_genomicGC12[record.id][1] = gc12_total_count

                    virus_genomicCpG[record.id][0] = cpg_count
                    virus_genomicCpG[record.id][1] = cpg_total_count
                    virus_genomicGpC[record.id][0] = gpc_count
                    virus_genomicGpC[record.id][1] = gpc_total_count
                    virus_genomicGpG[record.id][0] = gpg_count
                    virus_genomicGpG[record.id][1] = gpg_total_count
                    virus_genomicCpC[record.id][0] = cpc_count
                    virus_genomicCpC[record.id][1] = cpc_total_count
                else:
                    virus_genomicGC[record.id][0] += gc_count
                    virus_genomicGC[record.id][1] += gc_total_count
                    virus_genomicGC3[record.id][0] += gc3_count
                    virus_genomicGC3[record.id][1] += gc3_total_count
                    virus_genomicGC12[record.id][0] += gc12_count
                    virus_genomicGC12[record.id][1] += gc12_total_count

                    virus_genomicCpG[record.id][0] += cpg_count
                    virus_genomicCpG[record.id][1] += cpg_total_count
                    virus_genomicGpC[record.id][0] += gpc_count
                    virus_genomicGpC[record.id][1] += gpc_total_count
                    virus_genomicGpG[record.id][0] += gpg_count
                    virus_genomicGpG[record.id][1] += gpg_total_count
                    virus_genomicCpC[record.id][0] += cpc_count
                    virus_genomicCpC[record.id][1] += cpc_total_count

    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found. Please check the file path.")
    except Exception as e:
        print(f"An error occurred while processing {fasta_file}: {e}")

#print out genomic GC content
print(f"virus,genomincGC,genomicGC3,genomicGC12,genomicCpG,genomicGpC,genomicGpG,genomicCpC")
for virus_id, genomicGC in virus_genomicGC.items():
        gc=genomicGC[0]/genomicGC[1]*100
        gc3=virus_genomicGC3[virus_id][0]/virus_genomicGC3[virus_id][1]*100
        gc12=virus_genomicGC12[virus_id][0]/virus_genomicGC12[virus_id][1]*100
        gc_repeat=(virus_genomicGC3[virus_id][0]+virus_genomicGC12[virus_id][0])/(virus_genomicGC3[virus_id][1]+virus_genomicGC12[virus_id][1])*100
        cpg=virus_genomicCpG[virus_id][0]/virus_genomicCpG[virus_id][1]*1000
        gpc=virus_genomicGpC[virus_id][0]/virus_genomicGpC[virus_id][1]*1000
        gpg=virus_genomicGpG[virus_id][0]/virus_genomicGpG[virus_id][1]*1000
        cpc=virus_genomicCpC[virus_id][0]/virus_genomicCpC[virus_id][1]*1000
        print(f"{virus_id},{gc},{gc3},{gc12},{gc_repeat},{cpg},{gpc},{gpg},{cpc}")

# Write GC content and GC dinucleotide frequency in LIBSVM format
output_file = prefix_path + "gc_content_dinucleotide_libsvm.txt"
with open(output_file, "w") as out:
    for virus_id, segment_data in virus_segment_data.items():
        # Use a dummy label (e.g., negative -1) as we don't have actual labels
        label = -1
        # Construct the feature vector in LIBSVM format
        feature_list = []
        feature_index = 1
        for segment_index, data in enumerate(segment_data):
            if data is not None:
                gc_content, cpg_frequency, gpc_frequency, gpg_frequency, cpc_frequency = data
                feature_list.append(f"{feature_index}:{gc_content:.4f}")
                feature_index+=1
                feature_list.append(f"{feature_index}:{cpg_frequency:.4f}")
                feature_index+=1
                feature_list.append(f"{feature_index}:{gpc_frequency:.4f}")
                feature_index+=1
                feature_list.append(f"{feature_index}:{gpg_frequency:.4f}")
                feature_index+=1
                feature_list.append(f"{feature_index}:{cpc_frequency:.4f}")
                feature_index+=1
        # Write the line in LIBSVM format
        feature_vector = " ".join(feature_list)
        out.write(f"{label} {feature_vector}\n")

print(f"GC content and dinucleotide frequencies in LIBSVM format have been written to '{output_file}'.")
