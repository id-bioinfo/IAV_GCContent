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
        return 0  # Avoid division by zero
    return gc_count / total_count * 100

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
    return dinucleotide_count / total_count * 1000

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

prefix_path = "MinkH5_for_demo_using_cds/"
# Process each FASTA file
for i, fasta_file in enumerate(fasta_files, start=0):
    try:
        # Read the FASTA file
        with open(prefix_path + fasta_file, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                # Ensure the sequence is in uppercase
                sequence = str(record.seq).upper()  
                # Calculate GC content for this sequence
                gc_content = calculate_gc_content(sequence)
                # Calculate CpG frequency for this sequence
                cpg_frequency = calculate_dinucleotide_frequency(sequence, 'C', 'G')
                # Calculate GpC frequency for this sequence
                gpc_frequency = calculate_dinucleotide_frequency(sequence, 'G', 'C')
                # Calculate GpG frequency for this sequence
                gpg_frequency = calculate_dinucleotide_frequency(sequence, 'G', 'G')
                # Calculate CpC frequency for this sequence
                cpc_frequency = calculate_dinucleotide_frequency(sequence, 'C', 'C')

                # Store the data for this virus and segment (index i)
                virus_segment_data[record.id][i] = (gc_content, cpg_frequency, gpc_frequency, gpg_frequency, cpc_frequency)
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found. Please check the file path.")
    except Exception as e:
        print(f"An error occurred while processing {fasta_file}: {e}")

# Write GC content and GC dinucleotide frequency in LIBSVM format
output_file = "gc_content_dinucleotide_libsvm.txt"
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