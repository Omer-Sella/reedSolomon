# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:18:16 2024

@author: Omer
"""
import os, sys
reedSolomonProjectDir = os.environ.get('REEDSOLOMON')
if reedSolomonProjectDir == None: 
     reedSolomonProjectDir = "c:/users/omer/reedSolomon/reedSolomon/"
sys.path.insert(0, reedSolomonProjectDir)
from arithmetic import binaryFieldElement as galoisElement
from arithmetic import polynomial as polynomialClass
from arithmetic import generateExponentAndLogTables, gf128, generateExponentAndLogTables
from bchDecoder import bchDecoder
import numpy as np


def test_bchDecoder():
    encodedZeroData = np.zeros(126)
    eD, _ =  generateExponentAndLogTables()
    #print("Done generating log and exp dictionaries.")
    correctedVector, correctionVector, errorLocatorX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
    assert (np.all(correctedVector == 0))
    
def test_bchDecoder_single_bit_flip():
    encodedZeroData = np.zeros(126)
    eD, _ =  generateExponentAndLogTables()    
    for i in range(len(encodedZeroData)):
        encodedZeroData[i] = 1
        correctedVector, correctionVector, errorLocatorX = bchDecoder( receivedBinaryVecotor = encodedZeroData, 
                                                                      exponentDictionary = eD, 
                                                                      numberOfPowers = 16, 
                                                                      codewordLengthMaximal = 127)
        #errorLocatorX.printValues()
        #print(np.where(correctionVector == 1))
        # Error found in location i
        assert correctionVector[i] == 1 
        # Only one error was found
        assert (np.sum(correctionVector) == 1)
        # All errors were fixed
        assert (np.sum(correctedVector) == 0)
        # Reset the test vector
        encodedZeroData[i] = 0

def test_bchDecoder_single_bit_flip_order_check():
    encodedZeroData = np.zeros(126)
    eD, _ =  generateExponentAndLogTables()
    for i in range(len(encodedZeroData)):
        encodedZeroData[i] = 1
        correctedVector, correctionVector, errorLocatorX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
        assert errorLocatorX.order() == 1
        encodedZeroData[i] = 0
        
def coverage_bchDecoder_two_bit_flips():   
    """
    This function checks exhausitively, that all combinations of two bit flips are correctable.
    I removed it from being a test function because it takes too long.
    """
    from itertools import combinations
    from bchDecoder import syndromeCalculator
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData = np.zeros(126)
    #testCombinations = [ (7, 51), (7, 52), (7, 53), (7, 54), (7, 55), (7, 56), (7, 100),  (7, 101), (7, 102),]
    testCombinations = list(combinations(range(126),2))
    for pair in testCombinations:
        encodedZeroData[pair[0]] = 1
        encodedZeroData[pair[1]] = 1
        # Notice that the decoder needs to produce the error locator polynomial eX for this coverage !
        correctedVector, correctionVector, eX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
        receivedBinaryAsPolynomial = polynomialClass(coefficients = list(map(gf128, encodedZeroData)))
        syndromes = syndromeCalculator(eD, numberOfPowers = 16, receivedBinaryX = receivedBinaryAsPolynomial)
        # Error found in both locations
        assert correctionVector[pair[0]] == 1 
        assert correctionVector[pair[1]] == 1 
        # Only two errors were found
        assert (np.sum(correctionVector) == 2)
        # All errors were fixed
        assert (np.sum(correctedVector) == 0)
        print(pair)
        
        
        encodedZeroData[pair[0]] = 0
        encodedZeroData[pair[1]] = 0
        
def coverage_retrace_bug_error_in_first_coordinate(index):
    encodedZeroData = np.zeros(126)
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData[index] = 1
    # Notice that the decoder needs to produce the error locator polynomial eX for this coverage !
    correctedVector, correctionVector, eX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
    for i in eD.keys():
        print(eX.at(gf128(eD[i])).getValue())
    return correctedVector, correctionVector, eX



        

def test_connection_polynomial_for_two_errors_explicit_calculation_gf128():   
    """
    Explicit calculation using Todd K. Moon page 252
    """
    from itertools import combinations
    from bchDecoder import syndromeCalculator
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData = np.zeros(126)
    
    testCombinations = [ (7, 51), (7, 52), (7, 53), (7, 54), (7, 55), (7, 56), (7, 100),  (7, 101), (7, 102),]
    for pair in testCombinations:
        encodedZeroData[pair[0]] = 1
        encodedZeroData[pair[1]] = 1
        # Notice that the decoder needs to produce the error locator polynomial eX for this coverage !
        correctedVector, correctionVector, eX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
        receivedBinaryAsPolynomial = polynomialClass(coefficients = list(map(gf128, encodedZeroData)))
        syndromes = syndromeCalculator(eD, numberOfPowers = 16, receivedBinaryX = receivedBinaryAsPolynomial)
        gfOne = gf128(1)
        lambda1 = syndromes[0]
        lambda2 = (syndromes[2] + (syndromes[0] * syndromes[0] * syndromes[0])) / syndromes[0]
        explicitConnectionX = polynomialClass(coefficients = [lambda2, lambda1, gfOne])
        assert explicitConnectionX == eX
        encodedZeroData[pair[0]] = 0
        encodedZeroData[pair[1]] = 0
        
def test_connection_polynomial_for_three_errors_explicit_calculation_gf128():   
    """
    Explicit calculation using Todd K. Moon page 252
    """
    from itertools import combinations
    from bchDecoder import syndromeCalculator
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData = np.zeros(126)
    #testCombinations = np.random.choice(list(combinations(range(126),3)), 20)
    testCombinations = [ (0,7, 51), (1,7, 52), (2,7, 53), (3,7, 54), (4,7, 55), (5,7, 56), (6,7, 100),  (10,7, 101), (11,7, 102),]
    for errorLocations in testCombinations:
        for e in errorLocations:
            encodedZeroData[e] = 1
        # Notice that the decoder needs to produce the error locator polynomial eX for this coverage !
        correctedVector, correctionVector, eX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
        receivedBinaryAsPolynomial = polynomialClass(coefficients = list(map(gf128, encodedZeroData)))
        syndromes = syndromeCalculator(eD, numberOfPowers = 16, receivedBinaryX = receivedBinaryAsPolynomial)
        gfOne = gf128(1)
        lambda1 = syndromes[0]
        lambda2 = ((syndromes[1] * syndromes[1]) 
                   * syndromes[2] 
                   + syndromes[4]) / ((syndromes[0] * 
                                       syndromes[0] * 
                                       syndromes[0]) + 
                                       syndromes[2])
        lambda3 = ((syndromes[0] * 
                            syndromes[0] * 
                            syndromes[0]) + 
                            syndromes[2]) + syndromes[2] * lambda2
        explicitConnectionX = polynomialClass(coefficients = [lambda3, lambda2, lambda1, gfOne])
        if eX != explicitConnectionX:
            eX.printValues()
            explicitConnectionX.printValues()
        assert explicitConnectionX == eX
        for e in errorLocations:
            encodedZeroData[e] = 0
        
        
def test_connection_polynomial_for_four_errors_explicit_calculation_gf128():   
    """
    Explicit calculation using Todd K. Moon page 252
    """
    from itertools import combinations
    from bchDecoder import syndromeCalculator
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData = np.zeros(126)
    #testCombinations = np.random.choice(list(combinations(range(126),4)), 20)
    testCombinations = [ (0,7, 51,110), (1,7, 52,111), (2,7, 53,112), (3,7, 54,120), (4,7, 55,121), (5,7, 56,122), (6,7, 100,123),  (10,7, 101,124), (11,7, 102,125)]
    for errorLocations in testCombinations:
        for e in errorLocations:
            encodedZeroData[e] = 1
        # Notice that the decoder needs to produce the error locator polynomial eX for this coverage !
        correctedVector, correctionVector, eX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
        receivedBinaryAsPolynomial = polynomialClass(coefficients = list(map(gf128, encodedZeroData)))
        s = syndromeCalculator(eD, numberOfPowers = 16, receivedBinaryX = receivedBinaryAsPolynomial)
        gfOne = gf128(1)
        lambda1 = s[0]
        lambda2 = (s[0] * (s[6] + s[0]*s[0]*s[0]*s[0]*s[0]*s[0]*s[0]) + s[2] * (s[0]*s[0]*s[0]*s[0]*s[0] + s[4])) / (s[2] * (s[0]*s[0]*s[0] + s[2]) + s[0] * (s[0]*s[0]*s[0]*s[0]*s[0] + s[4]))
        lambda3 = s[0]*s[0]*s[0] + s[2] + s[0] * lambda2
        lambda4 =  (s[4] + s[0]*s[0]*s[2] + lambda2 * (s[0]*s[0]*s[0] + s[2])) / s[0]
        explicitConnectionX = polynomialClass(coefficients = [lambda4, lambda3, lambda2, lambda1, gfOne])
        assert explicitConnectionX == eX
        for e in errorLocations:
            encodedZeroData[e] = 0
            
if __name__ == "__main__":
    test_bchDecoder_single_bit_flip()