import gzip
from Bio import SeqIO
import os
import numpy as np
import matplotlib.pyplot as plt

def read_fastq(file_path):
    read_counts = {}
    with gzip.open(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)
            if seq in read_counts:
                read_counts[seq] += 1
            else:
                read_counts[seq] = 1
    return read_counts

def process_files(file_list):
    all_counts = {}
    for file_path in file_list:
        counts = read_fastq(file_path)
        file_name = os.path.basename(file_path)
        all_counts[file_name] = counts
    return all_counts

def calculate_tc_scaling_factors(all_counts):
    total_counts = {file_name: sum(counts.values()) for file_name, counts in all_counts.items()}
    average_total_count = sum(total_counts.values()) / len(total_counts)
    scaling_factors = {file_name: total / average_total_count for file_name, total in total_counts.items()}
    return scaling_factors

def calculate_uq_scaling_factors(all_counts):
    upper_quartiles = {}
    for file_name, counts in all_counts.items():
        nonzero_counts = [count for count in counts.values() if count > 0]
        upper_quartile = np.percentile(nonzero_counts, 75)
        upper_quartiles[file_name] = upper_quartile
    average_upper_quartile = sum(upper_quartiles.values()) / len(upper_quartiles)
    scaling_factors = {file_name: uq / average_upper_quartile for file_name, uq in upper_quartiles.items()}
    return scaling_factors

def normalize_counts(all_counts, scaling_factors):
    normalized_counts = {}
    for file_name, counts in all_counts.items():
        normalized_counts[file_name] = {read: count / scaling_factors[file_name] for read, count in counts.items()}
    return normalized_counts

def calculate_median_scaling_factors(all_counts):
    medians = {}
    for file_name, counts in all_counts.items():
        nonzero_counts = [count for count in counts.values() if count > 0]
        median = np.median(nonzero_counts)
        medians[file_name] = median
    average_median = sum(medians.values()) / len(medians)
    scaling_factors = {file_name: median / average_median for file_name, median in medians.items()}
    return scaling_factors

def plot_read_distribution(counts, title):
    for file_name, read_counts in counts.items():
        values = list(read_counts.values())
        plt.figure(figsize=(10, 6))
        plt.hist(values, bins=50, alpha=0.7)
        plt.title(f"{title} - {file_name}")
        plt.xlabel("Read Count")
        plt.ylabel("Frequency")
        plt.yscale('log')
        plt.show()