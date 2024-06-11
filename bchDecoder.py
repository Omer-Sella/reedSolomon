# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:13:47 2024

@author: Omer
"""

# A BCH decoder using key equation solver, chien search
from keyEquationSolver import keyEquationSolver
from arithmetic import binaryFieldElement as galoisElement
from arithmetic import polynomial as polynomialClass
from arithmetic import gf128
from arithmetic import generateExponentAndLogTables
import numpy as np
import time

def syndromeCalculator(exponentDictionary, numberOfPowers, receivedBinaryX):
    syndromes = []
    for i in range(numberOfPowers):
        # exponent table starts with gf(1) at index 0 (!), but s0 is receivedBinaryX(\alpha) !
        newSyndrome = receivedBinaryX.at(gf128(exponentDictionary[1 + i])) 
        syndromes.append(newSyndrome)
    return syndromes



def bchDecoder(receivedBinaryVecotor, exponentDictionary, numberOfPowers, codewordLengthMaximal):
    correctionVector = np.zeros(codewordLengthMaximal, dtype = np.int32)
    receivedBinaryX = polynomialClass(coefficients = list(map(gf128, receivedBinaryVecotor))) #list(map(gf128, receivedBinaryVecotor[::-1])))
    #receivedBinaryX.printValues()
    # Calculate syndromes
    syndromes = syndromeCalculator(exponentDictionary, numberOfPowers, receivedBinaryX)
    # Calculate error locator polynomial
    errorLocatorX = keyEquationSolver(polynomialClass, gf128, syndromes)
    #errorLocatorX.printValues()
    # Chien search 
    for i in range(0, len(exponentDictionary.keys())):#codewordLengthMaximal):
        # If \alpha ^i is a root of the error locator polynomial, then log( (\alpha ^i )^-1 ) is a location of an error.
        if (errorLocatorX.at(gf128(exponentDictionary[i]))) == 0: #Note: if we always use gf128, then cast the entire dictionary in advance
            # Since in GF(2^k) the multiplicative group is of order (2^k)-1, then the above value is just (codewordLengthMaximal - i) BUT (!!!) in IEEE notation the polynomials are "leading coefficient is leftmost" so it's just i ...
            #correctionVector[codewordLengthMaximal - i] = 1
            correctionVector[i] = 1
    correctedVector = (receivedBinaryVecotor + correctionVector[0 : len(receivedBinaryVecotor)]) %2
    #return correctedVector, correctionVector
    return correctedVector, correctionVector, errorLocatorX #For debug, communicate the error locator as well
            
            
        
        