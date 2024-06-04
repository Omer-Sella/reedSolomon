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


def bchDecoder(receivedBinaryVecotor, exponentDictionary, numberOfPowers, codewordLengthActual, codewordLengthMaximal):
    correctionVector = np.zeros(codewordLengthActual, dtype = np.int32)
    receivedBinaryX = polynomialClass(coefficients = list(receivedBinaryVecotor))
    
    # Calculate syndromes
    syndromes = []
    for i in range(numberOfPowers):
        #print("exponent dictionary at location" + str(i) + " is "+ str(exponentDictionary[i]))
        helper = gf128(exponentDictionary[i])
        print("The class of the eval point is:")
        print(helper.__class__)
        start =        time.time()
        
        newSyndrome = receivedBinaryX.at(helper)
        end = time.time()
        print("Time it took to evaluate at point is: " + str(end-start))
        print("Syndrome class is")
        print(newSyndrome.__class__)
        syndromes.append(newSyndrome)
        #print(syndromes)
    # Calculate error locator polynomial
    errorLocatorX = keyEquationSolver(polynomialClass, gf128, syndromes)
    # Chien search 
    for i in range(codewordLengthActual):
        if (errorLocatorX.at(exponentDictionary[i])).isZro():
            correctionVector[codewordLengthMaximal - i] = 1
    correctedVector = receivedBinaryVecotor + correctionVector
    return correctedVector, correctionVector
            
            
        
        