import numpy as np
import textwrap
import gzip
import random

class RefSeq(dict):
    """
    Reference sequence/s object - subclass of dict
    Format = {header:DNA(seq)}
    """

    def load_ref_file(self, in_file):
        """
        Load a FASTA-formatted reference file.  RNA is converted to DNA.
        :param in_file: path to FASTA formated reference file (compressed or uncompressed)
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
                    if header:
                        self[header] = DNA("".join(sequence_lines)).u_to_t()
                        sequence_lines = []
                    header = line[1:]
                else:
                    sequence_lines.append(line)
            if header:
                self[header] = DNA("".join(sequence_lines)).u_to_t()

    def seq_stats(self):
        seq_lens = []
        for seq in self.values():
            seq_lens.append(len(seq))
        return len(seq_lens), sum(seq_lens), np.mean(seq_lens)

    def kmer(self, kmer_len):
        for i in self.values():
            for j in range(len(i)-kmer_len+1):
                yield(DNA(i[j:j+kmer_len]))
    
    def write_ref_file(self, out_file):
        with open(out_file, 'w') as f:
            for header, seq in self.items():
                f.write(">{}\n".format(header))
                for line in textwrap.wrap(seq, 120):
                    f.write("{}\n".format(line))

    def gc_content(self):
        """
        Calculate and return the combined GC content of all the loaded sequences.
        The GC content is represented as a value between 0 and 1 and rounded to two decimal places.
        """
        total_bases = 0
        gc_bases = 0
        for seq in self.values():
            total_bases += len(seq)
            gc_bases += seq.count('G') + seq.count('C')

        if total_bases == 0:
            return 0.0  # Return 0.0 if there are no bases to avoid division by zero.
        return round(gc_bases / total_bases, 2)
    
    
tab = str.maketrans("ACTG", "TGAC")

class DNA(str):
    """
    DNA object
    """

    def __new__(cls, content):
        return super(DNA, cls).__new__(cls, content.upper())

    def __getitem__(self, item):
        return DNA(super(DNA, self).__getitem__(item))

    def __add__(self, item):
        return DNA(super(DNA, self).__add__(item))

    def split(self, sep=None, maxsplit=-1):
        return [DNA(i) for i in str.split(sep, maxsplit)]

    def u_to_t(self):
        return DNA(self.replace("U", "T"))

    def reverse_complement(self):
        """
        Reverse complement a DNA string
        :return: string of ACGTN
        """
        return self.translate(tab)[::-1]
    
    def canonicalise(self):
        a = self.reverse_complement()
        if a < self: 
            return a
        else:
            return self
        
class RandomSeqGen:
    @staticmethod
    def nonspecific_dsrna(length, gc_content=0.4):
        half_length = length / 2
        gcs = gc_content * half_length
        ats = half_length - gcs
        if not gcs.is_integer() or not ats.is_integer():
            raise ValueError(
                "It's impossible to generate a sequence with exactly the given length and GC content."
            )
        gcs = int(gcs)
        ats = int(ats)
        seq = list(gcs * "G" + gcs * "C" + ats * "A" + ats * "T")
        random.shuffle(seq)
        rand_seq = "".join(seq)
        return rand_seq