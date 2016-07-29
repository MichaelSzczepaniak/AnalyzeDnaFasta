

def readFasta(inFilePath=".\dna.example.fasta") :
    """ Reads a FASTA file inFilePath into a dictionary where
    the sequence identifiers are the keys and the DNA
    sequence for that record are the values.
    
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
        if line.startswith('>') :  # header line? Could also use line =='>'
            words = line.split()   # Split on space
            name = words[0][1:]    # Use everything right of > as dict key
            dnaSeqs[name] = ''     # Add new key to dict (assign value later)
        else :             # Line is not header, so append to the dna
                           # sequence string until we get to end of sequence.
            dnaSeqs[name] = dnaSeqs[name] + line
    f.close()
    
    return dnaSeqs