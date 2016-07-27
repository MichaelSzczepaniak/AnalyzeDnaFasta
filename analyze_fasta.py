#!/usr/bin/python3
import sys, getopt
import read_fasta as rf, dna_nrepeats as nreps

def usage() :
    print """
    analyze_fasta.py : reads a FASTA file and builds a
    dictionary with all sequences bigger than a given
    length.
    
    After reading the file a number of analysis can be run
    on the FASTA data.  These analysis include:
    
    1) tbd
    2) tbd
    
    analyze_fasta.py [-h] [--record-count] [--longest-seq] [shortest-seq]
                     [] [] []
    
    
    <filename>
    
    -h          print this message
    
    --record-count  output the number of FASTA records in <filename>
    
    
    
    <filename>      FASTA data file to be read and processed
    
    """

def main() :
    """
    Example of how to execute:
    From command line:  python dnaseqfasta.py dna2.fasta 6
    From ipython shell: run dnaseqfasta.py dna2.fasta 6
    
    
    
    """
    file_fasta = sys.argv[1]     # 1st arg should be the FASTA file
    n_repeat = int(sys.argv[2])  # 2nd arg should be the length of repeats
    data_fasta = rf.readFasta(file_fasta)