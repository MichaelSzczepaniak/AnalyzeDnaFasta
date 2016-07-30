#!/usr/bin/python3


def getNRepeats(dna, n=2, values_are_counts=True) :
    """ Returns a dictionary where:
    keys = lower case substrings of length n in dna and values are are either
    1) the frequency/count of each substring if values_are_counts = True
    (default) or 2) the index of each substring found in dna if
    values_are_counts = False.
    
    dna - string, rep'n of DNA sequence of A's, T's, G's and C's which
          can be upper or lower case
    n - integer, allowable value: 1 to len(dna), length of substrings
        to count frequencies on
    """
    dna = dna.lower()
    # Get all the substrings and set them as keys in our dictionary.
    repeats = {}
    for i in range(0, len(dna)-n) :
        substr = dna[i:i+n]
        if values_are_counts :
            repeats[substr] = 0
        else :
            repeats[substr] = []
    # Count or get start index of each instance of each key.
    for key in repeats.keys() :
        for j in range(0, len(dna)-n) :
            subdna = dna[j:j+n]
            if key == subdna :
                if values_are_counts :
                    repeats[key] += 1
                else :
                    repeats[key].append(j+1) # convert to 1-based indices
                
    return repeats
    
def getFirstMostFrequentRepeatN(dna_seqs, nrep) :
    """ Returns a string that is the first instance of the most frequently
    occuring repeat of size nrep.
    
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

def getAllMostFrequentRepeatN(dna_seqs, nreps) :
    """ Returns a 2-tuple where:
    tuple[0] = list of strings that are all the instances of the most
               frequently occuring repeats of size nrep.
    tuple[1] = integer number of occurances of the repeats in tuple[0]
    
    dna_seqs - Dictionary where keys = sequence identifiers and values =
    DNA seqeunce for that record.
    nrep - Size of the repeat to search on.
    """
    longest_n_repeats = []
    longest_n_sub_count = -1
    for seq_id in dna_seqs.keys() :
        dna = dna_seqs[seq_id]
        repeats = getNRepeats(dna, nreps)  # Get all repeats size n.
        for nrepeat, nrep_count in repeats.items() :
            if nrep_count > longest_n_sub_count :
                longest_n_sub_count = nrep_count
                longest_n_repeats = [nrepeat, ]
            elif nrep_count == longest_n_sub_count :
                longest_n_repeats.append(nrepeat)
    
    return (longest_n_repeats, longest_n_sub_count)
    
def getOccsOfRepeatInSingleSeq(dna, repeat) :
    """ Returns an integer which is the number of occurance of repeat
    in dna
    """
    repeat = repeat.lower()
    all_repeats = getNRepeats(dna, len(repeat))
    occurances = all_repeats.get(repeat, 0)
    
    return occurances
    
def getOccsOfRepeatInAllSeqs(dna_seqs, repeat) :
    """ Returns a dictionary with keys that are sequence ids and
    values that are integers which are the number of occurance of repeat
    in each of the sequences in dna_seqs
    
    dna_seqs - Dictionary where keys = sequence identifiers and values =
    DNA seqeunce for that record.
    repeat - string of A's, T's, G's and C's which representing a nucleotide
             sequence
    """
    result = {}
    for sid, dna in dna_seqs.items() :
        result[sid] = getOccsOfRepeatInSingleSeq(dna, repeat)
        
    return result
    
def getAllOccsOfRepeat(dna_seqs, repeat) :
    """ Returns an integer that's the count of all instance of repeat
    over all sequences in dna_seqs.
    
    dna_seqs - Dictionary where keys = sequence identifiers and values =
    DNA seqeunce for that record.
    repeat - string of A's, T's, G's and C's which representing a nucleotide
             sequence
    """
    result = 0
    allOccsOfRepeat = getOccsOfRepeatInAllSeqs(dna_seqs, repeat)
    for count in allOccsOfRepeat.values() :
        result += count
    
    return result
    
