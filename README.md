# DnaSeqFastA

## Description
This is small Python project that reads a file containing DNA sequences in a multi-FASTA format and does the following analysis:

1. Determines the number of records in the file.
2. Determines the longest and shortest sequences and creates a histogram of sequence lengths
3. Determines the frequencies of the longest and shortest sequences and provide their sequence identifiers.
4. Determines the longest open reading frame (ORF) in the FASTA file
5. Determines the sequence identifier containing the longest ORF in the FASTA file
6. Determines the longest ORF for a given sequence identifier
7. Determines the starting position of the longest ORF for a given sequence identifier
8. Determine how many times a given repeat occurs in a given FASTA file
9. Determine the most frequent repeat in a given set of repeats.

## Usage
To run the program, from the teminal: <pre>python dnaseqfasta.py fastafile repeatsequencefile outfile</pre>

### Inputs
The first parameter, **fastafile** is a FASTA file which will be parsed and analyzed for the quantities described above.
The second parameter, **repeatsequencefile** is a file containing a single line of 'A', 'T', 'G', and 'C' characters representing nucleotide sequences.  These characters may be lower case (as permitted by [NCBI in doing BLAST queries](http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml).

### Outputs
The program outputs a file call **dnaFastaAnalysis.txt** to the directory from which **dnaseqfasta.py** is executed.

## Version
This project was developed using the Anaconda distribution of Python 3.5.1.
