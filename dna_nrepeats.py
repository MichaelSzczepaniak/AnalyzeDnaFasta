#!/usr/bin/python3
import sys, read_fasta as rfa

def getNRepeats(dna, n=2, values_as_locations=True) :
    """ Returns a dictionary where the keys are the substrings
    of length n in dna and values are are either 1) the index of each
    substring found in dna if values_as_locations = True or
    2) the frequency/count of each substring if values_as_locations = False.
    
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
        if values_as_locations :
            repeats[substr] = []
        else :
            repeats[substr] = 0
    # Count the instances of each keys
    for key in repeats.keys() :
        for j in range(0, len(dna)-n) :
            subdna = dna[j:j+n]
            if key == subdna :
                if values_as_locations :
                    repeats[key].append(j+1) # convert to 1-based indices
                else :
                    repeats[key] += 1
                
    return repeats
    
def getOccsOfRepeatInSingleSeq(dna_seq, repeat) :
    """ Returns an integer which is the number of occurance of repeat
    in dna_seq
    """
    all_repeats = getNRepeats(dna, len(repeat), False)
    occurances = all_repeats.get(repeat, 0)
    
    return occurances
    
def getOccsOfRepeatInAllSeqs(dna_seqs, repeat) :
    """ Returns a dictionary with keys that are the keys of dna_seqs and
    values that are integers which are the number of occurance of repeat
    in each of the sequences in dna_seqs
    """
    result = {}
    for sid, dna in dna_seqs.items() :
        result[sid] = getMaxOccOfRepeatInSingleSeq(dna, repeat)
        
    return result
    
def getAllOccsOfRepeat(dna_seqs, repeat) :
    """ Returns an integer that's the count of all instance of repeat
    in each sequence of dna_seqs.
    """
    result = 0
    repeats_in_all_seqs = getOccsOfRepeatInAllSeqs(dna_seqs, repeat)
    for count in repeats_in_all_seqs.values() :
        result += count
    
    return result
    
