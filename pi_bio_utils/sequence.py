import numpy as np
import textwrap
import gzip
import random
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns

class RefSeq(dict):
    """
    A subclass of the Python dict class for storing reference DNA sequences.

    Attributes:
        dict: The parent class for storing reference sequences.

    The format for storing reference sequences in the dictionary is {header: DNA(seq)}
    """

    def load_ref_file(self, in_file):
        """
        Load a FASTA-formatted reference file. RNA is converted to DNA.

        Args:
            in_file (str): Path to FASTA formatted reference file (compressed or uncompressed).
        """
        if in_file.endswith(".gz"):
            file_opener = gzip.open
        else:
            file_opener = open

        with file_opener(in_file, "rt") as f:
            header = ""
            sequence_lines = []
            for line in f:
                line = line.strip()
                if line == "":
                    pass
                elif line[0] == ">":
                    if header and sequence_lines:
                        seq = DNA("".join(sequence_lines)).u_to_t()
                        if seq:  # Check if the sequence is not empty.
                            self[header] = seq
                        sequence_lines = []
                    header = line[1:]
                else:
                    sequence_lines.append(line)
            if header and sequence_lines:
                seq = DNA("".join(sequence_lines)).u_to_t()
                if seq:  # Check if the sequence is not empty.
                    self[header] = seq

    def load_multiple_alignment(
        self, in_file, start=None, end=None, truncate_header=False
    ):
        """
        Load a multiple sequence alignment file and slice sequences between
        start and end positions, while removing dashes.

        Args:
            in_file (str): Path to FASTA-formatted multiple sequence alignment file.
            start (int, optional): The starting position to slice the sequence. 1-indexed.
            end (int, optional): The ending position to slice the sequence. 1-indexed.
            truncate_header (bool, optional): If True, truncate the header at the first space.

        Note:
            If start and end are provided, the slice is inclusive.
        """
        if in_file.endswith(".gz"):
            file_opener = gzip.open
        else:
            file_opener = open

        with file_opener(in_file, "rt") as f:
            header = ""
            sequence_lines = []
            for line in f:
                line = line.strip()
                if line == "":
                    pass
                elif line[0] == ">":
                    if header and sequence_lines:
                        seq = DNA("".join(sequence_lines))
                        if start is not None and end is not None:
                            seq = seq[start - 1 : end]
                        elif start is not None:
                            seq = seq[start - 1 :]
                        elif end is not None:
                            seq = seq[:end]
                        seq = DNA(seq.replace("-", ""))
                        if seq:
                            self[header] = seq
                        sequence_lines = []
                    header = line[1:]
                    if truncate_header:
                        header = header.split(" ")[0]
                else:
                    sequence_lines.append(line)

            # Add the final sequence
            if header and sequence_lines:
                seq = DNA("".join(sequence_lines))
                if start is not None and end is not None:
                    seq = seq[start - 1 : end]
                elif start is not None:
                    seq = seq[start - 1 :]
                elif end is not None:
                    seq = seq[:end]
                seq = DNA(seq.replace("-", ""))
                if seq:
                    self[header] = seq

    def seq_stats(self):
        """
        Calculate and return the number, sum of lengths, and mean length of the loaded reference sequences.

        Returns:
            tuple: A tuple containing the number of sequences, sum of lengths, and mean length.
        """
        seq_lens = []
        for seq in self.values():
            seq_lens.append(len(seq))
        return len(seq_lens), sum(seq_lens), np.mean(seq_lens)

    def kmer(self, kmer_len):
        """
        Generate k-mers of specified length from the loaded reference sequences.

        Args:
            kmer_len (int): Length of k-mers to generate.

        Yields:
            DNA: k-mer of the specified length.
        """
        for i in self.values():
            for j in range(len(i) - kmer_len + 1):
                yield (DNA(i[j : j + kmer_len]))

    def write_ref_file(self, out_file):
        """
        Write the loaded reference sequences to a FASTA-formatted file.

        Args:
            out_file (str): Path to the output FASTA-formatted file.

        Note:
            Each sequence line will be wrapped at 120 characters.
        """
        with open(out_file, "w") as f:
            for header, seq in self.items():
                f.write(">{}\n".format(header))
                for line in textwrap.wrap(seq, 120):
                    f.write("{}\n".format(line))

    def gc_content(self):
        """
        Calculate and return the combined GC content of all the loaded sequences.

        Returns:
            float: The GC content represented as a value between 0 and 1, rounded to two decimal places.

        Note:
            If there are no bases, the function will return 0.0 to avoid division by zero.
        """
        total_bases = 0
        gc_bases = 0
        for seq in self.values():
            total_bases += len(seq)
            gc_bases += seq.count("G") + seq.count("C")

        if total_bases == 0:
            return 0.0  # Return 0.0 if there are no bases to avoid division by zero.
        return round(gc_bases / total_bases, 2)


tab = str.maketrans("ACTG", "TGAC")


class DNA(str):
    """
    A class for DNA strings based on Python's built-in string class.

    Attributes:
        str: The parent class for storing DNA sequences.

    The DNA class inherits all string methods but returns DNA objects instead of strings where applicable.
    The class ensures the DNA sequence is uppercase and provides additional methods specific to DNA sequences.
    """

    def __new__(cls, content):
        """
        Create a new DNA instance with uppercase content.

        Args:
            content (str): The DNA sequence to store in the object.

        Returns:
            DNA: A new DNA object with uppercase content.
        """
        return super(DNA, cls).__new__(cls, content.upper())

    def __getitem__(self, item):
        """
        Retrieve an item or slice from the DNA object.

        Args:
            item (int, slice): The index or slice to retrieve.

        Returns:
            DNA: The item or slice as a DNA object.
        """
        return DNA(super(DNA, self).__getitem__(item))

    def __add__(self, item):
        """
        Concatenate a string or DNA object to the DNA object.

        Args:
            item (str, DNA): The string or DNA object to concatenate.

        Returns:
            DNA: A new DNA object containing the concatenated sequence.
        """
        return DNA(super(DNA, self).__add__(item))

    def split(self, sep=None, maxsplit=-1):
        """
        Split the DNA sequence by a separator.

        Args:
            sep (str, optional): The separator to split by. If None, split by whitespace. Defaults to None.
            maxsplit (int, optional): The maximum number of splits. Defaults to -1, which means no limit.

        Returns:
            list: A list of DNA objects representing the split sequence.
        """
        return [DNA(i) for i in super(DNA, self).split(sep, maxsplit)]

    def u_to_t(self):
        """
        Replace all instances of uracil (U) with thymine (T) in the DNA sequence.

        Returns:
            DNA: A new DNA object with uracil replaced by thymine.
        """
        return DNA(self.replace("U", "T"))

    def reverse_complement(self):
        """
        Return the reverse complement of the DNA sequence.

        Returns:
            DNA: A new DNA object representing the reverse complement of the sequence.

        Note:
            The reverse complement is formed by reversing the DNA sequence and replacing each base with its complement.
        """
        return self.translate(tab)[::-1]

    def canonicalise(self):
        """
        Return the lexicographically smaller of the DNA sequence and its reverse complement.

        Returns:
            DNA: The lexicographically smaller DNA object between the original sequence and its reverse complement.
        """
        a = self.reverse_complement()
        if a < self:
            return a
        else:
            return self


class RandomSeqGen:
    """
    A class for generating random DNA sequences.

    This class provides static methods for generating DNA sequences with specified characteristics, such as GC content.
    """

    @staticmethod
    def nonspecific_dsrna(length, gc_content=0.4):
        """
        Generate a non-specific double-stranded RNA (dsRNA) sequence with a specified length and GC content.

        Args:
            length (int): The length of the dsRNA sequence to generate.
            gc_content (float, optional): The desired GC content as a value between 0 and 1. Defaults to 0.4.

        Returns:
            str: The generated dsRNA sequence.

        Raises:
            ValueError: Raised if it's impossible to generate a sequence with exactly the given length and GC content.

        Note:
            The function creates a dsRNA sequence with the given GC content by calculating the number of G, C, A, and T bases.
            It then shuffles the sequence randomly.
        """
        half_length = length / 2
        gcs = gc_content * half_length
        ats = half_length - gcs
        if not gcs.is_integer() or not ats.is_integer():
            raise ValueError(
                "It's impossible to generate a sequence with exactly the given length and GC content providing count(C) == count(G)."
            )
        gcs = int(gcs)
        ats = int(ats)
        seq = list(gcs * "G" + gcs * "C" + ats * "A" + ats * "T")
        random.shuffle(seq)
        rand_seq = "".join(seq)
        return rand_seq

<<<<<<< HEAD
=======
class ReadLengthDistribution:
    def __init__(self, treatment_to_files, min_len, max_len):
        self.treatment_to_files = treatment_to_files
        self.min_len = min_len
        self.max_len = max_len
        self.all_counts = defaultdict(list)
        
        # Calculate counts automatically when a new instance is created
        self.calculate_counts()

    def calculate_counts(self):
        for treatment, fastq_gz_paths in self.treatment_to_files.items():
            for fastq_gz_path in fastq_gz_paths:
                count_by_length = Counter()
                total_reads_in_range = 0

                with gzip.open(fastq_gz_path, 'rt') as f:
                    while True:
                        identifier = f.readline().strip()
                        sequence = f.readline().strip()
                        plus = f.readline().strip()
                        quality = f.readline().strip()

                        if not identifier:
                            break

                        read_length = len(sequence)

                        if self.min_len <= read_length <= self.max_len:
                            total_reads_in_range += 1
                            count_by_length[read_length] += 1

                if total_reads_in_range > 0:
                    for read_length in count_by_length:
                        count_by_length[read_length] = (count_by_length[read_length] / total_reads_in_range) * 1e6

                self.all_counts[treatment].append(count_by_length)

    def plot(self):
        plt.figure(figsize=(12, 8))

        color_palette = sns.color_palette("bright", len(self.all_counts))
        sns.set_palette(color_palette)

        marker_styles = ['o', 's', 'v', '^', '<', '>', 'p', 'P', '*', 'h', 'H', 'D', 'd', 'X']

        all_lengths = set()

        for idx, (treatment, list_of_counts) in enumerate(self.all_counts.items()):
            aggregate_counts = defaultdict(list)

            for counts in list_of_counts:
                for read_length, count in counts.items():
                    aggregate_counts[read_length].append(count)
                    all_lengths.add(read_length)

            mean_counts = {}
            sem_counts = {}
            for read_length, counts in aggregate_counts.items():
                mean_counts[read_length] = np.mean(counts)
                sem_counts[read_length] = np.std(counts) / np.sqrt(len(counts))

            sorted_lengths = sorted(mean_counts.keys())
            sorted_means = [mean_counts[length] for length in sorted_lengths]
            sorted_sems = [sem_counts[length] for length in sorted_lengths]

            marker = marker_styles[idx % len(marker_styles)]

            plt.errorbar(sorted_lengths, sorted_means, yerr=sorted_sems,
                         label=treatment, color=color_palette[idx], marker=marker)

        plt.xticks(sorted(all_lengths))
        plt.xlabel("Read Length")
        plt.ylabel("Reads Per Million (RPM)")
        plt.title("Read Length Distribution by Treatment")
        plt.legend()
        plt.show()


>>>>>>> fdb68b3ebe373a1f4bc1846112977fb8e3222668
class BigWigHelper:
    def __init__(self, in_csv):
        """
        Initialize the BigWigHelper class with the input file path and read the data.
        
        Parameters:
        - in_csv: str, path to the input CSV file. Expected format: "chromosome,length"
        """
        self.in_csv = in_csv
        self.results = {}
        
        try:
            with open(self.in_csv, 'r') as f_in:
                next(f_in)  # skip header
                for line in f_in:
                    parts = line.split(",")
                    chrom = parts[0]
                    length = int(parts[1])  # convert to integer, assuming length is an integer
                    self.results.setdefault(chrom, length)
        except FileNotFoundError:
            print(f"File {self.in_csv} not found.")
        except Exception as e:
            print(f"An error occurred: {e}")
    
    def write_sorted(self, out_txt):
        """
        Sorts the chromosome names and their lengths, and writes the sorted data to a TXT file.
        
        Parameters:
        - out_txt: str, path to the output TXT file.
        
        Returns:
        None
        """
        sorted_dict = dict(sorted(self.results.items()))
        
        try:
            with open(out_txt, 'w') as f_out:
                for k, v in sorted_dict.items():
                    f_out.write(f"{k}\t{v}\n")
        except FileNotFoundError:
            print(f"File {out_txt} not found.")
        except Exception as e:
            print(f"An error occurred: {e}")