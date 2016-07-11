import sys

def main():
    print('main function')

if __name__ == "__main__": main()  # allow func calls before def's

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
        else :
            print('not stop codon: ', codon)
            
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
    except IOError:
        print("File doesn't exist!  Exiting.")
        
        
    dnaSeqs = {}    # init empty dict
    for line in f:  # iterate thru lines in file
        line = line.rstrip()       # discard newlines
        if line.startswith('>') :  # Are we on a header line?
            words = line.split()   # split on space
            name = words[0][1:]    # Use everything right of > as dict key
            dnaSeqs[name] = ''        # Add new key to dict (assigne value later)
        else :             # Line is not header, seq append DNA sequence
                           # Cont. appending parts of the seq
            dnaSeqs[name] = dnaSeqs[name] + line
    f.close()
    
    return dnaSeqs