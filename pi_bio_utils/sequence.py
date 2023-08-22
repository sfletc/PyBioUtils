import numpy as np
import textwrap
import gzip
import random


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
        if in_file.endswith('.gz'):
            file_opener = gzip.open
        else:
            file_opener = open

        with file_opener(in_file, 'rt') as f:
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
