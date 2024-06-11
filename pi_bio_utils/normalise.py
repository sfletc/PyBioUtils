import gzip
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write
import os
import sys

def parse_fastq(file_path):
    sequences = {}
    with gzip.open(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)
            if seq in sequences:
                sequences[seq] += 1
            else:
                sequences[seq] = 1
    return sequences

# def parse_fastq_files(file_paths):
#     all_sequences = []
#     for file_path in file_paths:
#         sequences = parse_fastq(file_path)
#         for seq, count in sequences.items():
#             all_sequences.append({'sequence': seq, 'count': count, 'file': os.path.basename(file_path)})
#     return pd.DataFrame(all_sequences)

def normalize_rpm(sequences_df):
    total_reads = sequences_df['count'].sum()
    sequences_df['normalized'] = sequences_df['count'] / total_reads * 1e6
    return sequences_df

def normalize_median(sequences_df):
    median_count = sequences_df['count'].median()
    sequences_df['normalized'] = sequences_df['count'] / median_count
    return sequences_df

def normalize_quantile(sequences_df):
    ranks = sequences_df['count'].rank(method='min')
    sorted_counts = sequences_df['count'].sort_values().values
    quantile_normalized_counts = [sorted_counts[int(rank - 1)] for rank in ranks]
    sequences_df['normalized'] = quantile_normalized_counts
    return sequences_df

def normalize_upper_quartile(sequences_df):
    upper_quartile = sequences_df['count'].quantile(0.75)
    sequences_df['normalized'] = sequences_df['count'] / upper_quartile
    return sequences_df

def write_fasta(file_path, sequences_df):
    sequences_df = sequences_df.sort_values(by='normalized', ascending=False).reset_index(drop=True)
    records = []
    for index, row in sequences_df.iterrows():
        header = f"read_{index+1}-{row['normalized']:.2f}"
        sequence = row['sequence']
        record = SeqRecord(Seq(sequence), id=header, description="")
        records.append(record)
    write(records, file_path, "fasta")

def plot_distribution(data, title, output_dir):
    plt.figure(figsize=(10, 6))

    # Plot before normalization
    for file_name, data in data.items():
        sorted_counts = data['count'].sort_values()
        plt.plot(sorted_counts.values, label=f'Raw counts {file_name}')
    
    plt.xlabel('Read Rank')
    plt.ylabel('Read Count')
    plt.yscale('log')
    plt.title(title)
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.savefig(os.path.join(output_dir, "before_normalization.png"), bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()


def process_file(file_path, all_sequences):
    sequences = parse_fastq(file_path)
    for seq, count in sequences.items():
        all_sequences.append({'sequence': seq, 'count': count, 'file': os.path.basename(file_path)})
    return sequences

def apply_rpm_normalization(file_paths, output_dir):
    before_data = {}
    after_data = {}
    for file_path in file_paths:
        file_name = os.path.basename(file_path).replace('.fastq.gz', '').replace('.fq.gz', '')
        sequences_df = pd.DataFrame(parse_fastq(file_path).items(), columns=['sequence', 'count'])
        before_data[file_name] = sequences_df.copy()
        sequences_df = normalize_rpm(sequences_df)
        after_data[file_name] = sequences_df
        output_file = os.path.join(output_dir, f"{file_name}_normalized.fasta")
        write_fasta(output_file, sequences_df)
    return before_data, after_data

def apply_combined_normalization(file_paths, sequences_df, normalization_method, output_dir):
    before_data = {}
    after_data = {}

    # Non-RPM normalization
    if normalization_method == "median":
        sequences_df = normalize_median(sequences_df)
    elif normalization_method == "quantile":
        sequences_df = normalize_quantile(sequences_df)
    elif normalization_method == "upper_quartile":
        sequences_df = normalize_upper_quartile(sequences_df)
    else:
        raise ValueError("Unsupported normalization method")

    # Split and write to individual FASTA files
    for file_path in file_paths:
        file_name = os.path.basename(file_path).replace('.fastq.gz', '').replace('.fq.gz', '')
        file_sequences_df = sequences_df[sequences_df['file'] == file_name].copy()
        before_data[file_name] = file_sequences_df.copy()
        after_data[file_name] = file_sequences_df
        output_file = os.path.join(output_dir, f"{file_name}_normalized.fasta")
        write_fasta(output_file, file_sequences_df)

    return before_data, after_data

def main(input_files, normalization_method, output_dir):
    # Split the comma-separated list of input files
    file_paths = input_files.split(',')
    
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    if normalization_method == "RPM":
        normalised = apply_rpm_normalization(file_paths, output_dir)
    else:
        # Data for combined normalization
        all_sequences = []
        for file_path in file_paths:
            process_file(file_path, all_sequences)
        
        # Create dataframe from combined data
        sequences_df = pd.DataFrame(all_sequences)
        
        normalised = apply_combined_normalization(file_paths, sequences_df, normalization_method, output_dir)

    # Plot distributions
    plot_distribution(normalised, output_dir)