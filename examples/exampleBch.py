# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 09:27:50 2024

@author: Omer
"""
import os, sys
reedSolomonProjectDir = os.environ.get('REEDSOLOMON')
if reedSolomonProjectDir == None: 
     reedSolomonProjectDir = "c:/users/omer/reedSolomon/reedSolomon/"
sys.path.insert(0, reedSolomonProjectDir)
import numpy as np
from arithmetic import generateExponentAndLogTables, polynomial, gf128, binaryFieldElement as gf2
from bchDecoder import bchDecoder

def example_binarySymmetricChannel():
    # Create the generator polynomial with coefficients [1, 0, 1, 0, 0, 1,    1,    1,0 ,1, 0,1, 0, 1, 0, 1, 1] as a polynomial over gf128
    gX = polynomial( coefficients = list(map(gf2, np.array([1, 0, 1, 0, 0, 1,    1,    1,0 ,1, 0,1, 0, 1, 0, 1, 1]) )))
    #
    binaryData = np.random.randint(0,2,110)
    # cast the data into gf128 and make into polynomial over gf128
    data = polynomial(coefficients = list(map(gf2, binaryData)))
    # Calculate the parity as the remainder modulu gX, of the lift of the data by the order of gX
    data = data.lift(gX.order())
    parity = data.modulu(gX)
    data.coefficients[len(data.coefficients) - gX.order() : ] = parity.coefficients
    encodedBinaryData = np.array([data.coefficients[i].value for i in range(len(data.coefficients))])
    # Flip 2 bits
    errorLocations = np.random.choice(126,2,replace = False)
    encodedBinaryData[errorLocations] = 1 - encodedBinaryData[errorLocations]
    
    eD, _ =  generateExponentAndLogTables()
    correctedVector, correctionVector, errorLocatorX = bchDecoder( receivedBinaryVecotor = encodedBinaryData,
                                                                  exponentDictionary = eD,
                                                                  numberOfPowers = 16,
                                                                  codewordLengthMaximal = 127)
    errorLocatorX.printValues()
    
    
def example_cacheABinaryEncoder():
    """
    generate a binary matrix G, 
    which, as a linear transformation over GF(2), 
    is identical to encoding using polynomial division with the generator polynomial 
    # g(x) =      x16 +  x14 +    x11 + x10 + x9 + x7 + x5 +  x3 +  x+ 1
    gX = [1, 0, 1, 0, 0, 1,    1,    1,0 ,1, 0,1, 0, 1, 0, 1, 1]
    """
    # Create the generator polynomial with coefficients [1, 0, 1, 0, 0, 1,    1,    1,0 ,1, 0,1, 0, 1, 0, 1, 1] as a polynomial over gf128
    gX = polynomial( coefficients = list(map(gf2, np.array([1, 0, 1, 0, 0, 1,    1,    1,0 ,1, 0,1, 0, 1, 0, 1, 1]) )))
    onez = np.eye(110)
    zeroPadding = np.zeros((16,110))
    parityMatrix = np.vstack((onez, zeroPadding))
    for i in range(110):
        data = polynomial(coefficients = list(map(gf2, parityMatrix[:,i])))
        parity = data.modulu(gX)
        parityBits = np.array([parity.coefficients[i].value for i in range(len(parity.coefficients))])
        print(parityBits)
        parityMatrix[110 : 126, i] = parityBits
    return parityMatrix
    
    
    
                      