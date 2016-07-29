
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
        #print("  i=", i, ", codon=", codon, "|")
        if codon in check_codons :
            #print("  ", codon_type, "CODON found!\n")
            next_start_index = i
            break
        #else :
            #print("  ", codon_type, "codon NOT found...\n")
    
    return next_start_index

def getOpenReadingFrames(raw_dna, reading_frame=1) :
    """ Returns a list called orfs_list containing 2-tuples. Each tuple
    corresponds to an Open Reading Frame (ORF) where:
    tuple[0] = index of start codon
    tuple[1] = length of the ORF
    
    raw_dna - string, rep'n of DNA sequence of A's, T's, G's and C's
          which can be upper or lower case
    reading_frame - int, valid values: 1, 2, or 3 corresponding to reading
          frames 1, 2, or 3 respespectively to be used in reading dna string
    
    orfs_list is sorted by tuple[1] (ORF length) ascending order
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
            true_start_codon_index = start_codon_index+reading_frame-1
            orf_record = (true_start_codon_index, orf_length)
            orfs_list.append(orf_record)
            #print("  *********** adding tuple:", orf_record, " ***********\n")
        else :
            break  # No more stop codons following start codons: exit
        # Look for next start codon.
        start_codon_python_index = \
        getNextCodonIndex(dna, start_codon_python_index + 3)
    # Sort list by third element in tuple: (how does lambda thingy work?...)
    # http://stackoverflow.com/questions/3121979/how-to-sort-list-tuple-of-lists-tuples
    orfs_list = sorted(orfs_list, key=lambda tup: tup[1])
    
    return orfs_list
    
def getShortLongestORFsInAll(dna_dict, want_shortest=True, reading_frame=1) :
    """ Returns a 3-tuple where 
    
    tuple[0] - length of the ORF, shortest if want_shortest = True or
               longest if want_shortest = False
    tuple[1] - index of the start codon for the ORF of interest.
    tuple[2] - seq_id where the ORF was found.
    
    If no ORFs could be found, return tuple will have -1 in first 2 elements.
    
    dna - string, rep'n of DNA sequence of A's, T's, G's and C's which
          can be upper or lower case
    reading_frame - reading frame to find longest or shortest ORF on. Valid
                    values are 1, 2, 3, or 0 if over all reading frames
    """
    minmaxlen = -1
    minmaxindex = -1
    minmaxseq = "unassigned"
    if not(reading_frame in (0, 1, 2, 3)) : return minmaxlen
    all_rfs = (1, 2, 3)  # All valid reading frames
    rframes = [reading_frame, ] if (reading_frame in all_rfs) else all_rfs
    print("rframes = {}".format(rframes))
    for rframe in rframes :
        for seq_id, dna in dna_dict.items() :
            orfs = getOpenReadingFrames(dna, rframe)
            if len(orfs) < 1 :
                # print("len={}, index={}, seq={}".format(orfs[-1][1], orfs[-1][0], seq_id))
                continue  # Couldn't find an ORF in this seq_id
            if want_shortest :
                if orfs[0][1] < minmaxlen :
                    minmaxindex = orfs[0][0]
                    minmaxlen = orfs[0][1]
                    minmaxseq = seq_id
            else :
                if orfs[-1][1] > minmaxlen :
                    minmaxindex = orfs[-1][0]
                    minmaxlen = orfs[-1][1]
                    minmaxseq = seq_id
        
    return (minmaxlen, minmaxindex, minmaxseq)
    
def getLengthLongestORF(dna_dict, seq_id, reading_frame=1) :
    """ Returns a 2-tuples:
    tuple[0] - length of the longest ORF in seq_id
    tuple[1] - reading for which the longeset ORF is computed
    
    dna_dict - keys = sequence ids, values = DNA sequence of A's, T's, G's
               and C's which can be upper or lower case
    """
    dna = dna_dict[seq_id]
    if reading_frame in (1, 2, 3) :
        orfs = getOpenReadingFrames(dna, reading_frame)
        if len(orfs) > 1 :
            return (orfs[-1][0], reading_frame)
        else :
            return (0, reading_frame)
    elif reading_frame == 0 :
        lorf = -1  # Init length of longest ORF
        for rf in (1, 2, 3) :
            orfs = getOpenReadingFrames(dna, rf)
            if orfs[-1][1] > lorf :
                lorf = orfs[-1][1]
                rframe = rf
        return(lorf, rf)
    else :
        print("getLengthLongestORF expecting reading_frame to be 1,2,3, or 0.")
        return (-1, reading_frame)
    
