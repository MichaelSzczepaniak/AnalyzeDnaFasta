#!/usr/bin/python3
import sys, read_fasta as rfa

def getNRepeatCounts(dna, n=2) :
    """ Returns a dictionary where the keys are the substrings
    of length n in dna and values are the frequency/count of each
    substring found in dna.
    
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
        repeats[substr] = 0
    # Count the instances of each keys
    for key in repeats.keys() :
        for j in range(0, len(dna)-n) :
            subdna = dna[j:j+n]
            if key == subdna :
                repeats[key] += 1
                
    return repeats
    
def getNRepeatPositions(dna, n=2) :
    """ Returns a dictionary where the keys are the substrings
    of length n in dna and values are the index of each substring
    found in dna.
    
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
        repeats[substr] = []
    # Count the instances of each keys
    for key in repeats.keys() :
        for j in range(0, len(dna)-n) :
            subdna = dna[j:j+n]
            if key == subdna :
                repeats[key].append(j)
                
    return repeats