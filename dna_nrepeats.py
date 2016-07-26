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
                    repeats[key].append(j)
                else :
                    repeats[key] += 1
                
    return repeats