# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:13:47 2024

@author: Omer
"""

# A BCH decoder using key equation solver, chien search
from keyEquationSolver import keyEquationSolver
from arithmetic import binaryFieldElement as galoisElement
from arithmetic import polynomial
from arithmetic import generateExponentAndLogTables
import numpy as np
import time

def syndromeCalculator(exponentDictionary, gfType, numberOfPowers, receivedBinaryX):
    syndromes = []
    for i in range(numberOfPowers):
        # exponent table starts with gf(1) at index 0 (!), but s0 is receivedBinaryX(\alpha) ! 
        newSyndrome = receivedBinaryX.at(gfType(exponentDictionary[i + 1]))
        syndromes.append(newSyndrome)
    return syndromes

def forneyCalculator(syndromeX, errorLocatorX, t):
    """
    Input:
        syndromeX - syndrome polynomial, i.e., polynomial which coefficients are the syndromes.
        errorLocatorX - polynomial which was found by the key equation solver
        t - half the designed distance of the code
        
    Output:
        omegaX - polynomial that satisfies omegaX = syndromeX * errorLocatorX modulo x^(2t)
    """
    omegaX = (syndromeX * errorLocatorX).ignoreFromDegree(2*t-1) #Remember that ignoreFromDegree(d) ignores everything with degree strictly greater than d
    return omegaX


def bchDecoder(receivedBinaryVecotor, gfType, exponentDictionary, numberOfPowers, codewordLengthMaximal, reedSolomon = False):
    correctionVector = np.zeros(codewordLengthMaximal, dtype = np.int32)
    receivedBinaryX = polynomial(coefficients = list(map(gfType, receivedBinaryVecotor))) #list(map(gf128, receivedBinaryVecotor[::-1])))
    # Calculate syndromes
    syndromes = syndromeCalculator(exponentDictionary, gfType, numberOfPowers, receivedBinaryX)
    # Calculate error locator polynomial
    errorLocatorX = keyEquationSolver(polynomial, gfType, syndromes)
    #errorLocatorX.printValues()
    # Chien search 
    # Note the problem here: exponentDictionary starts at index 0, and in index 0 it says 1 (not alpha)
    # So to begin chien search at \alpha we start at ecponentDictionary[1]. We check all the way to len(exponentDictionary.keys() + 1) so, say 127, 
    # but when we want to access the exponentDictionary at 127 we take the 127 modulus to actually check at 0.
    for i in range(1, len(exponentDictionary) + 1): 
        # If \alpha ^i is a root of the error locator polynomial, then log( (\alpha ^i )^-1 ) is a location of an error.
        if (errorLocatorX.at(gfType(exponentDictionary[i % (len(exponentDictionary))]))) == 0: #Note: there is a potential speed up here: cast the entire dictionary in advance to the appropriate data type (say gf128)
            # Since in GF(2^k) the multiplicative group is of order (2^k)-1, then the above value is just (codewordLengthMaximal - i) BUT (!!!) in IEEE notation the polynomials are "leading coefficient is leftmost" so it's just i ...
            #correctionVector[codewordLengthMaximal - i] = 1
            if reedSolomon== False:
                correctionVector[i - 1] = 1
            else:
                omegaX = forneyCalculator(polynomial(syndromes), errorLocatorX, numberOfPowers / 2)
                correctionVector[i - 1] = (omegaX.at(exponentDictionary[i % (len(exponentDictionary))])) / (errorLocatorX.d().at(exponentDictionary[i % (len(exponentDictionary))]))
    correctedVector = (receivedBinaryVecotor + correctionVector[0 : len(receivedBinaryVecotor)]) %2
    return correctedVector, correctionVector, errorLocatorX #For debug, communicate the error locator as well
            
            
        
        