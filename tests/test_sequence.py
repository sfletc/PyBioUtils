import unittest
from pi_bio_utils.sequence import RefSeq, DNA, RandomSeqGen

class TestRefSeq(unittest.TestCase):

    def setUp(self):
        self.refseq = RefSeq()
        self.refseq.load_ref_file("tests/test.fa")

    def test_load_ref_file(self):
        self.assertEqual(len(self.refseq), 2)
        self.assertEqual(self.refseq["1"], "ATCACTACTGCA")
        self.assertEqual(self.refseq["3"], "ATCGATCGCAATCGATCGAT")

    def test_seq_stats(self):
        num_seqs, total_length, mean_length = self.refseq.seq_stats()
        self.assertEqual(num_seqs, 2)
        self.assertEqual(total_length, 32)
        self.assertAlmostEqual(mean_length, 16, 1)

    def test_kmer(self):
        kmer_len = 4
        kmers = list(self.refseq.kmer(kmer_len))
        self.assertEqual(len(kmers), 26)
        self.assertEqual(kmers[0], "ATCA")
        self.assertEqual(kmers[-1], "CGAT")

    def test_gc_content(self):
        gc_content = self.refseq.gc_content()
        self.assertAlmostEqual(gc_content, 0.44, 2)

    def test_write_ref_file(self):
        self.refseq.write_ref_file("tests/test_output.fa")
        with open("tests/test_output.fa", "r") as f:
            lines = f.readlines()
            self.assertEqual(lines[0].strip(), ">1")
            self.assertEqual(lines[1].strip(), "ATCACTACTGCA")
            self.assertEqual(lines[2].strip(), ">3")
            self.assertEqual(lines[3].strip(), "ATCGATCGCAATCGATCGAT")


class TestDNA(unittest.TestCase):

    def test_init(self):
        dna = DNA("acgt")
        self.assertEqual(str(dna), "ACGT")

    def test_getitem(self):
        dna = DNA("ACGT")
        self.assertEqual(dna[1], "C")
        self.assertEqual(dna[1:3], "CG")

    def test_add(self):
        dna1 = DNA("ACGT")
        dna2 = DNA("TGCA")
        dna3 = dna1 + dna2
        self.assertEqual(str(dna3), "ACGTTGCA")

    def test_split(self):
        dna = DNA("ACGTGCTG")
        split_dna = dna.split("G")
        self.assertEqual([str(d) for d in split_dna], ["AC", "T", "CT", ""])

    def test_u_to_t(self):
        dna = DNA("AUCGU")
        self.assertEqual(str(dna.u_to_t()), "ATCGT")

    def test_reverse_complement(self):
        dna = DNA("ACGT")
        self.assertEqual(str(dna.reverse_complement()), "ACGT")  # Assuming you have implemented the 'tab' variable for translation in the reverse_complement method.

    def test_canonicalise(self):
        dna = DNA("ACGT")
        rev_comp = DNA("ACGT")
        self.assertEqual(str(dna.canonicalise()), "ACGT")
        self.assertEqual(str(rev_comp.canonicalise()), "ACGT")

class TestRandomSeqGen(unittest.TestCase):

    def test_nonspecific_dsrna(self):
        # Test if the length of the sequence generated is correct.
        seq = RandomSeqGen.nonspecific_dsrna(12, 0.5)
        self.assertEqual(len(seq), 12)
        
        # Test if the GC content of the sequence generated is correct.
        gc_content = (seq.count("G") + seq.count("C")) / len(seq)
        self.assertAlmostEqual(gc_content, 0.5, delta=0.1)

        # Test if an error is raised when generating a sequence with impossible GC content.
        with self.assertRaises(ValueError):
            RandomSeqGen.nonspecific_dsrna(10, 0.75)

if __name__ == '__main__':
    unittest.main()

