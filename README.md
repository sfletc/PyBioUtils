Here's a basic `README.md` file template for your project:

```markdown
# PyBioUtils

PyBioUtils is a Python package that provides tools for DNA sequence manipulation. This package includes three main classes: `RefSeq`, `DNA`, and `RandomSeqGen`.

## Features

- Load a FASTA-formatted reference file.
- Calculate the statistics of a DNA sequence.
- Generate k-mers of a given length from a DNA sequence.
- Write a reference sequence to a FASTA file.
- Calculate the GC content of a DNA sequence.
- Convert uracil (U) to thymine (T) in a DNA sequence.
- Generate the reverse complement of a DNA sequence.
- Generate a non-specific double-stranded RNA sequence.

## Installation

1. Clone the repository:

```bash
git clone https://github.com/sfletc/PyBioUtils.git
cd PyBioUtils
```

2. Install the package:

```bash
pip install .
```

## Usage

### RefSeq

```python
from pi_bio_utils.sequence import RefSeq

refseq = RefSeq()
refseq.load_ref_file("path/to/reference.fasta")
stats = refseq.seq_stats()
gc_content = refseq.gc_content()
```

### DNA

```python
from pi_bio_utils.sequence import DNA

seq = DNA("ATCGATCGATCG")
rev_complement = seq.reverse_complement()
```

### RandomSeqGen

```python
from pi_bio_utils.sequence import RandomSeqGen

random_seq = RandomSeqGen.nonspecific_dsrna(length=100, gc_content=0.4)
```

## Dependencies

- numpy
- textwrap3
- gzip-reader (optional, for gzip support)

## Author

Stephen Fletcher

## License

This project is licensed under the MIT License.

## Contact

You can contact the author at steveofbrisbane@gmail.com.
```

