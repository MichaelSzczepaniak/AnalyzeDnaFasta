#!/usr/bin/python3
import sys
import numpy as np

def main():
    """
    dnaseqfasta.py dna.test01.fasta repseq.txt
    """
    file_fasta = sys.argv[1]   # 1st arg should be the FASTA file
    file_repseq = sys.argv[2]  # 2nd arg should be the seq to check for repeats
    data_fasta = readFasta(file_fasta)
    record_count = getRecordCount(data_fasta)
    input_line = "input file = {}".format(file_fasta)
    output_content = [input_line, ]
    output_content.append("record count = {}".format(record_count))
    longest_sequences = getShortLongSeqs(data_fasta, shortest=False)
    output_content.extend(getSeqLengthContent(longest_sequences, "longest"))
    shortest_sequences = getShortLongSeqs(data_fasta, shortest=True)
    output_content.extend(getSeqLengthContent(shortest_sequences, "shortest"))
    writeOutputFile(output_content)

def getSeqLengthContent(length_tuple, length_type="longest") :
    """
    """
    return_list = ["----------",
        (length_type + " sequence = {}".format(length_tuple[0])),
        (length_type + " sequence count = {}".format(length_tuple[1])),
        (length_type + " sequence ids:")]
    for seq_id in length_tuple[2] :
        return_list.append(seq_id)
        
    return(return_list)
    
def writeOutputFile(content_list) :
    outfile = open('dnaFastaAnalysis.txt', 'w')
    for line in content_list :
        print(line, file = outfile, end = '\n')
    outfile.close()

def getRecordCount(fasta_dat) :
    """ fasta_dat - dictionary with keys = sequence ids and
    values = string that represents a DNA sequence of A's T's, G's and C's
    """
    return len(fasta_dat.keys())
    
def getShortLongSeqs(fasta_data, shortest=True) :
    """ Returns a 3-tuple where the first elements is the length of the
    shortest or longest sequence, the second element is the number of 
    sequences that have the shortest or longest length, and the third element
    is a list of sequence ids of sequences that have this shortest or longest
    length.
    
    fasta_data - a dictionary with keys that are sequence ids and values that
                 are sequences
    shortest - boolean: If True, all outputs are with respect to the shortest
               length sequence(s). If False, outputs are wrt to longest
               length sequence(s).
    """
    id_seqs = []      # Store the ids of the shortest or longest seq's
    # Init length of shortest or longest sequence
    length_seq = sys.maxsize if shortest else -sys.maxsize
    for seq_id, dna_seq in sorted(fasta_data.items()) :
        seq_len = len(fasta_data[seq_id])  # Get length of seq for this id
        print("seq_id:", seq_id, "\nhas length=", seq_len)
        # If this sequence is the shortest or longest we've seen so far,
        # update length_seq and reinit id_seqs
        if(shortest) :
            if(seq_len < length_seq) :
                length_seq = seq_len
                id_seqs = [seq_id,]
                print('new shortest sequence =', seq_id, '\nhas length = ', seq_len)
            elif(seq_len == length_seq) :
                id_seqs.append(seq_id)
                print('adding new seq to shortest list:\n', seq_id)
        else :
            if(seq_len > length_seq) :
                length_seq = seq_len
                id_seqs = [seq_id,]
                print('new longest sequence =', seq_id, '\nhas length = ', seq_len)
            elif(seq_len == length_seq) :
                id_seqs.append(seq_id)
                print('adding new seq to longest list:\n', seq_id)
            
    return((length_seq, len(id_seqs), id_seqs))

# test string
# dna = "atgtaaatatgctagatgcccat"
#        01234567890123456789012
#           *        *
# expected list: [3, 12]
def get_stop_codon(dna, frame=0) :
    """
    Returns a list of integers that are the indices of stop codons present
    in the string dna.
    
    dna - string, rep'n of DNA sequence of A's, T's, G's and C's which
          can be lower case
    frame - int, valid values: 0, 1, or 2 corresponding to reading frames
            1, 2, or 3 respespectively to be used in reading dna string
    checks if given dna sequence contains a stop codon
    """
    stop_codons_found = [-1,]
    stop_codons=['tga', 'tag', 'taa'] # nomalize to lower case
    for i in range(frame, len(dna), 3) :
        codon = dna[i:i+3].lower()  # codon triplet to check
        if codon in stop_codons :
            if stop_codons_found[0] < 0 :
                stop_codons_found[0] = i
            else :
                stop_codons_found.append(i)
        #else :
        #    print('not stop codon: ', codon)
            
    return stop_codons_found


def readFasta(inFilePath=".\dna.example.fasta") :
    """ Reads a FASTA file inFilePath into a dictionary where
    the sequence identifiers are the keys and the DNA
    seqeunce for that record are the values.
    
    Keyword aguments:
    inFilePath -- path to the FASTA file to be read
                  (default dna.example.fasts in current dir)
    """
    try:
        f = open(inFilePath)
        print("FASTA file read successfully.")
    except IOError:
        print("File doesn't exist!  Exiting.")
        sys.exit(0)
        
    dnaSeqs = {}    # init empty dict
    for line in f:  # iterate thru lines in file
        line = line.rstrip()       # discard newlines
        if line.startswith('>') :  # Are we on a header line?
            words = line.split()   # Split on space
            name = words[0][1:]    # Use everything right of > as dict key
            dnaSeqs[name] = ''     # Add new key to dict (assigne value later)
        else :             # Line is not header, seq append DNA sequence
                           # Cont. appending parts of the seq
            dnaSeqs[name] = dnaSeqs[name] + line
    f.close()
    
    return dnaSeqs
    
if __name__ == "__main__": main()  # allow func calls before def's