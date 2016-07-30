# AnalyzeDnaFasta

## Description
This is a small Python project that reads a file containing DNA sequences in a multi-FASTA format and does the following analysis:

1. Determines the number of records in the file.
2. Determines the longest and shortest sequences.
3. Determines the frequencies of the longest and shortest sequences and provide their sequence identifiers.
4. Determines the longest open reading frame (ORF) in the FASTA file.
5. Determines the sequence identifier containing the longest ORF in the FASTA file
6. Determines the longest ORF for a given sequence identifier.
7. Determines the starting position of the longest ORF for a given sequence identifier.
8. Determine how many times a given repeat occurs in a given FASTA file.
9. Determine the most frequent repeat in a given set of repeats.

## Usage
To run the program, from the command line, execute: <pre>python analyze\_fasta.py <fastafile> </pre>

Example: <pre>python analyze\_fasta.py data/dna.example1.fasta --record\_count</pre>

### Inputs
The first parameter, **fastafile** is always to location of the FASTA file which will be parsed and analyzed for the quantities described above.
The second parameter, **option** is always one of the following options:  

<ul>
<li><pre>--record_count</pre> Return the number of records in the input FASTA file.</li>
<li><pre>--longest_seq</pre> Return the length of the longest sequence in the input FASTA file.</li>
<li><pre>--shortest_seq</pre> Return the length of the shortest sequence in the input FASTA file.</li>
<li><pre>--longest_orf n</pre> Return the longest open reading frame (ORF) for reading frame n = 1, 2, or 3. Set n = 0 to obtain the longest ORF of any reading frame.</li>
</ul>

These characters may be lower case (as permitted by [NCBI in doing BLAST queries](http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml)).

### Outputs
The program outputs the requested results to the terminal.

## Python Version
This project was developed using the Anaconda distribution of Python 3.5.1.
