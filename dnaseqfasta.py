#!/usr/bin/python3
import sys
import numpy as np

def main():
    """
    run dnaseqfasta.py dna.test01.fasta repseq.txt
    run dnaseqfasta.py dna.test02.fasta repseq.txt
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
    # Get the ORFs in each record of data
    output_content.append("----------")
    for seq_id in data_fasta.keys() :
        print("**********START processing sequence = ", seq_id, " **********")
        output_content.append("sequence id = {}".format(seq_id))
        dna_seq = data_fasta[seq_id]
        for rf_offset in [0, 1, 2] :
            reading_frame = rf_offset + 1
            print("   >>>>>>> READING FRAME = ", reading_frame, " <<<<<<<")
            orfs = getOpenReadingFrames(dna_seq, reading_frame)
            output_content.append("  number of ORF{}s = ".format(reading_frame))
            for orf in orfs :
                output_content.append("  length = {}, start = {}".format(orf[2], orf[0]))
            output_content.append("  ***")
        output_content.append("----------")
        print("**********END processing sequence = ", seq_id, " **********")
    
    writeOutputFile(output_content)

def getSeqLengthContent(length_tuple, length_type="longest") :
    """ Returns a list of strings which describe the results of
    computing the longest and shortest sequences.
    
    length_tuple - the 3 tuple returned from getShortLongSeqs
    """
    return_list = ["----------",
        (length_type + " sequence = {}".format(length_tuple[0])),
        (length_type + " sequence count = {}".format(length_tuple[1])),
        (length_type + " sequence ids:")]
    for seq_id in length_tuple[2] :
        return_list.append(seq_id)
        
    return(return_list)
    
def writeOutputFile(content_list, outputFileName='dnaFastaAnalysis.txt') :
    """ Writes item in content_list as a line to the file
    named outputFileName
    """
    outfile = open(outputFileName, 'w')
    for line in content_list :
        print(line, file = outfile, end = '\n')
    outfile.close()

def getRecordCount(fasta_dat) :
    """ fasta_dat - dictionary with keys = sequence ids and
    values = string that represents a DNA sequence of A's T's, G's and C's
    """
    return len(fasta_dat.keys())
    
def getShortLongSeqs(fasta_data, shortest=True) :
    """ Find the length of longest or shortest sequence.
    Returns a 3-tuple where:
    tuple[0] - length of the shortest or longest sequence
    tuple[1] - number of sequences that have the shortest or longest length
    tuple[2] - list of sequence ids of sequences that have this shortest or
               longest length.
    
    fasta_data - a dictionary with keys that are sequence ids and values that
                 are sequences
    shortest - boolean: If True, all outputs are with respect to the shortest
               length sequence(s). If False, outputs are wrt to longest
               length sequence(s).
    """
    id_seqs = []      # Store the ids of the shortest or longest seq's
    # Init length of shortest or longest sequence
    length_seq = sys.maxsize if shortest else -1
    for seq_id, dna_seq in sorted(fasta_data.items()) :
        seq_len = len(fasta_data[seq_id])  # Get length of seq for this id
        #print("seq_id:", seq_id, "\nhas length=", seq_len)
        # If this sequence is the shortest or longest we've seen so far,
        # update length_seq and reinit id_seqs
        if(shortest) :
            if(seq_len < length_seq) :
                length_seq = seq_len
                id_seqs = [seq_id,]
                #print('new shortest sequence =', seq_id, '\nhas length = ', seq_len)
            elif(seq_len == length_seq) :
                id_seqs.append(seq_id)
                #print('adding new seq to shortest list:\n', seq_id)
        else :
            if(seq_len > length_seq) :
                length_seq = seq_len
                id_seqs = [seq_id,]
                #print('new longest sequence =', seq_id, '\nhas length = ', seq_len)
            elif(seq_len == length_seq) :
                id_seqs.append(seq_id)
                #print('adding new seq to longest list:\n', seq_id)
            
    return((length_seq, len(id_seqs), id_seqs))

# test:
# dna = "ATGGGTATGGGGTGA", start codons at 0 and 6, stop codon at 12
# If start_index = 0, start codon at 0 found, or stop codon at 12 found.
# If start_index = 3, start codon at 6 found, or stop codon at 12 found.
# NO START OR STOP CODONS ARE FOUND if start_index is anything other
# than 0 or 3. THEREFORE, BE CAREFUL TO SET start_index properly!
def getNextCodonIndex(dna, start_index=0, find_start_codon=True,
                      check_codons=['tga', 'tag', 'taa']) :
    """ Returns the next index in dna starting from start_index that is a
    start codon (i.e. ATG or atg) if find_start_codon==True (default).
    
    dna - string, rep'n of DNA sequence of A's, T's, G's and C's which
          can be upper or lower case
    start_index - int, Python index to start sequence reads
    find_start_codon - boolean. If start_codon==False and nothing is passed
    in for check_codons, then the first index of a stop codon that is after
    start_index is returned.
    
    If start_codon==False and a check_codons list is passed in, then the
    first index of a codon that matches any check_codons element is returned.
    
    If no codon is found that meets the above conditions -1 is returned.
    """

    codon_type = "START"
    if not find_start_codon :
        if set(check_codons) == set(('TAA', 'TAG', 'TGA')) :
            codon_type = "STOP"
        elif set(check_codons) == set(('taa', 'tag', 'tga')) :
            codon_type = "STOP"
        else :
            codon_type = "CUSTOM"
    if find_start_codon : check_codons = ['atg', ] # Look for start codon index
    next_start_index = -1
    for i in range(start_index, len(dna), 3) :
        codon = dna[i:i+3].lower()
        print("  i=", i, ", codon=", codon, "|")
        if codon in check_codons :
            print("  ", codon_type, "CODON found!\n")
            next_start_index = i
            break
        else :
            print("  ", codon_type, "codon NOT found...\n")
    
    return next_start_index

def getOpenReadingFrames(raw_dna, reading_frame=1) :
    """ Returns a list called orfs_list containing 3-tuples. Each tuple
    corresponds to an Open Reading Frame (ORF) where:
    tuple[0] = index of start codon
    tuple[1] = index of stop codon
    tuple[2] = length of the ORF
    
    raw_dna - string, rep'n of DNA sequence of A's, T's, G's and C's
          which can be upper or lower case
    reading_frame - int, valid values: 1, 2, or 3 corresponding to reading
          frames 1, 2, or 3 respespectively to be used in reading dna string
    
    orfs_list is sorted by tuple[2] (ORF length) ascending order
    """
    dna = raw_dna[reading_frame-1:]  # Read only bases in the reading frame
    orfs_list = []
    start_codon_python_index = getNextCodonIndex(dna) # Start at beginning
    while start_codon_python_index > -1 :
        stop_codon_python_index = \
        getNextCodonIndex(dna, start_codon_python_index + 3, False)
        if stop_codon_python_index > -1 :
            orf_length = stop_codon_python_index + 3 - start_codon_python_index
            # Create returned tuples as 1-based indices.
            start_codon_index = start_codon_python_index + 1
            stop_codon_index = stop_codon_python_index + 1
            orf_record = (start_codon_index, stop_codon_index, orf_length)
            orfs_list.append(orf_record)
            print("  *********** adding tuple:", orf_record, " ***********\n")
        else :
            break  # No more stop codons following start codons: exit
        # Look for next start codon.
        start_codon_python_index = \
        getNextCodonIndex(dna, start_codon_python_index + 3)
    # Sort list by third element in tuple: (how does lambda thingy work?...)
    # http://stackoverflow.com/questions/3121979/how-to-sort-list-tuple-of-lists-tuples
    orfs_list = sorted(orfs_list, key=lambda tup: tup[2])
    
    return orfs_list

def getShortLongestORFs(dna, want_shortest=True, reading_frame=1) :
    """ Returns the length of the shortest Open Reading Frame ORF if
    want_shortest = True in dna. Otherwise, returns the longest ORF.
    dna - string, rep'n of DNA sequence of A's, T's, G's and C's which
          can be upper or lower case
    reading_frame - int, valid values: 1, 2, or 3 corresponding to reading
          frames 1, 2, or 3 respespectively to be used in reading dna string
    """
    orfs = getOpenReadingFrames(dna, reading_frame)
    if want_shortest :
        result = orfs[0][2]
    else :
        result = orfs[-1][2]
        
    return result
    
# test string
# dna = "atgtaaatatgctagatgcccat"
#        01234567890123456789012
#           *        *
# expected list: [3, 12]
def get_stop_codons(dna, reading_frame=1) :
    """
    Returns a list of integers that are the indices of stop codons present
    in the string dna.
    
    dna - string, rep'n of DNA sequence of A's, T's, G's and C's which
          can be upper or lower case
    reading_frame - int, valid values: 1, 2, or 3 corresponding to reading
          frames 1, 2, or 3 respespectively to be used in reading dna string
    
    Function checks if given dna sequence contains a stop codon
    """
    stop_codons_found = []
    stop_codons=['tga', 'tag', 'taa'] # nomalize to lower case
    reading_index = reading_frame - 1
    for i in range(reading_index, len(dna), 3) :
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
                  (default dna.example.fasta in current dir)
    """
    try:
        f = open(inFilePath)
        print("FASTA file read successfully.")
    except IOError:
        print("File doesn't exist!  Exiting.")
        sys.exit(0)
        
    dnaSeqs = {}    # init empty dict
    for line in f:  # iterate thru lines in file
        line = line.rstrip()       # remove trailing (right) white space (\n)
        if line.startswith('>') :  # Are we on a header line?
            words = line.split()   # Split on space
            name = words[0][1:]    # Use everything right of > as dict key
            dnaSeqs[name] = ''     # Add new key to dict (assign value later)
        else :             # Line is not header, so append to the dna
                           # sequence string until we get to end of sequence.
            dnaSeqs[name] = dnaSeqs[name] + line
    f.close()
    
    return dnaSeqs
    
if __name__ == "__main__": main()  # allow func calls before def's