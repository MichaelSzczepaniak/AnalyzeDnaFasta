#!/usr/bin/python3


def getNRepeats(dna, n=2, values_are_counts=True) :
    """ Returns a dictionary where the keys are the substrings
    of length n in dna and values are are either 1) the index of each
    substring found in dna if values_are_counts = False or
    2) the frequency/count of each substring if values_are_counts = True
    (default).
    
    dna - string, rep'n of DNA sequence of A's, T's, G's and C's which
          can be upper or lower case
    n - integer, allowable value: 1 to len(dna), length of substrings
        to count frequencies on
    """
    dna = dna.lower()
    # Get all the substrings and set them as keys in our dictionary
    repeats = {}
    for i in range(0, len(dna)-n) :
        substr = dna[i:i+n]
        if values_are_counts :
            repeats[substr] = 0
        else :
            repeats[substr] = []
    # Count the instances of each keys
    for key in repeats.keys() :
        for j in range(0, len(dna)-n) :
            subdna = dna[j:j+n]
            if key == subdna :
                if values_are_counts :
                    repeats[key] += 1
                else :
                    repeats[key].append(j+1) # convert to 1-based indices
                
    return repeats
    
def getMostFrequentRepeatN(dna_seqs, nrep) :
    """ Returns a string that is the most frequently occuring repeat
    of size nrep.
    
    dna_seqs - Dictionary where keys = sequence identifiers and values =
    DNA seqeunce for that record.
    nrep - Size of the repeat to search on.
    
    """
    longest_n_repeat = ""
    longest_n_sub_count = -1
    for seq_id in dna_seqs.keys() :
        dna = dna_seqs[seq_id]
        repeats = getNRepeats(dna, nrep)
        for nrepeat, nrep_count in repeats.items() :
            if nrep_count > longest_n_sub_count :
                longest_n_sub_count = nrep_count
                longest_n_repeat = nrepeat
    
    return longest_n_repeat
    
    
def getOccsOfRepeatInSingleSeq(dna, repeat) :
    """ Returns an integer which is the number of occurance of repeat
    in dna
    """
    all_repeats = getNRepeats(dna, len(repeat))
    occurances = all_repeats.get(repeat, 0)
    
    return occurances
    
def getOccsOfRepeatInAllSeqs(dna_seqs, repeat) :
    """ Returns a dictionary with keys that are sequence ids and
    values that are integers which are the number of occurance of repeat
    in each of the sequences in dna_seqs
    """
    result = {}
    for sid, dna in dna_seqs.items() :
        result[sid] = getOccsOfRepeatInSingleSeq(dna, repeat)
        
    return result
    
def getAllOccsOfRepeat(dna_seqs, repeat) :
    """ Returns an integer that's the count of all instance of repeat
    in each sequence in dna_seqs.
    """
    result = 0
    getAllOccsOfRepeat = getOccsOfRepeatInAllSeqs(dna_seqs, repeat)
    for count in getAllOccsOfRepeat.values() :
        result += count
    
    return result
    
