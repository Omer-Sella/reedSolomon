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
        assert correctionVector[i] == 1
        encodedZeroData[i] = 0

def test_bchDecoder_single_bit_flip_order_check():
    encodedZeroData = np.zeros(126)
    eD, _ =  generateExponentAndLogTables()
    for i in range(len(encodedZeroData)):
        encodedZeroData[i] = 1
        correctedVector, correctionVector, errorLocatorX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
        assert errorLocatorX.order() == 1
        encodedZeroData[i] = 0
        
def coverage_retrace_bug_error_in_first_coordinate(index):
    encodedZeroData = np.zeros(126)
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData[index] = 1
    # Notice that the decoder needs to produce the error locator polynomial eX for this coverage !
    correctedVector, correctionVector, eX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
    for i in eD.keys():
        print(eX.at(gf128(eD[i])).getValue())
    return correctedVector, correctionVector, eX
    

def coverage_example_8_8_tkmoon():
    pass

    
def test_connection_polynomial_for_two_errors_explicit_calculation_gf128():   
    """
    Explicit calculation using Todd K. Moon page 252
    """
    from itertools import combinations
    from bchDecoder import syndromeCalculator
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData = np.zeros(126)
    testCombinations = np.random.choice(combinations(range(126),2), 20)
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
    testCombinations = np.random.choice(combinations(range(126),3), 20)
    for errorLocations in testCombinations:
        encodedZeroData[errorLocations] = 1
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
        assert explicitConnectionX == eX
        encodedZeroData[errorLocations] = 0
        
def test_connection_polynomial_for_four_errors_explicit_calculation_gf128():   
    """
    Explicit calculation using Todd K. Moon page 252
    """
    from itertools import combinations
    from bchDecoder import syndromeCalculator
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData = np.zeros(126)
    testCombinations = np.random.choice(combinations(range(126),4), 20)
    for errorLocations in testCombinations:
        encodedZeroData[errorLocations] = 1
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
        encodedZeroData[errorLocations] = 0