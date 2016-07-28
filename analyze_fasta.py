#!/usr/bin/python3
import sys, argparse
import read_fasta as rf, dna_orfs as orf, dna_nrepeats as nreps

def main() :
    """
    Example of how to execute:
    From command line:  python analyze_fasta.py <fasta file> --option
    From ipython shell: run analyze_fasta.py <fasta file> --option
    where option can be:
    record_count - to return the number of records in the input FASTA file
    longest_seq - to return the length of the longest sequence in the input
                  FASTA file
    shortest_seq - to return the length of the shortest sequence in the input
                  FASTA file
    
    """
    file_fasta = sys.argv[1]     # FASTA file, required regardless of option
    
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', nargs=1)
    parser.add_argument("--record_count", action='store_true',
    help="Return the number of records in the input FASTA file")
    parser.add_argument("--longest_seq", action='store_true',
    help="Return the length of the longest sequence in the input FASTA file")
    parser.add_argument("--shortest_seq", action='store_true',
    help="Return the length of the longest sequence in the input FASTA file")
    
    parser.add_argument("--longest_orf", type=int, nargs=1,
    help="Return the longest open reading frame (ORF) for reading frames \
    1, 2, 3, or all (0)")
    
    args = parser.parse_args()
    if args.filename :
        data_fasta = rf.readFasta(file_fasta)  # Always read the FASTA file.
    else :
        print("FASTA data file needs to be specified")
        sys.exit(1)
    
    if args.record_count :
        print("Record count = {}".format(getRecordCount(data_fasta)))
    elif args.longest_seq :
        long_seq = getShortLongSeqs(data_fasta, shortest=False)
        print("                 Length of longest sequence = {}".\
        format(long_seq[0]))
        print("Count of sequences that have longest length = {}".\
        format(len(long_seq[2])))
        print("Sequences with longest length:")
        for seq in long_seq[2] :
            print(seq)
    elif args.shortest_seq :
        short_seq = getShortLongSeqs(data_fasta, shortest=True)
        print("                 Length of shortest sequence = {}".\
        format(short_seq[0]))
        print("Count of sequences that have shortest length = {}".\
        format(len(short_seq[2])))
        print("Sequences with longest length:")
        for seq in short_seq[2] :
            print(seq)
    elif args.longest_orf :
        rframe = args.longest_orf[0]
        lorf = orf.getShortLongestORFsInAll(data_fasta, False, rframe)
        if rframe in (1, 2, 3) :
            print("longest orf in reading frame {} = {}".format(rframe, lorf))
        elif rframe == 0 :
            print("longest orf in all reading frames = {}".format(lorf))
        else :
            print("Reading frame parameter must be 1, 2, 3, or 0!")
    else :
        print("Unimplemented option... TODO")

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
        # print("seq_id:", seq_id, "\nhas length=", seq_len)
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
                # print('new longest sequence =', seq_id, '\nhas length = ', seq_len)
            elif(seq_len == length_seq) :
                id_seqs.append(seq_id)
                # print('adding new seq to longest list:\n', seq_id)
            
    return((length_seq, len(id_seqs), id_seqs))
    
if __name__ == "__main__": main()  # Load all function in module, then run main()