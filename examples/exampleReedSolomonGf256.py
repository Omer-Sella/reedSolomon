# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 09:31:34 2024

@author: Omer Sella

### Reed Solomon code over GF(2^8) example
"""
import os, sys
reedSolomonProjectDir = os.environ.get('REEDSOLOMON')
if reedSolomonProjectDir == None: 
    #Note that this is an example, so the relevant modules will be in its parent directory, hence os.path.dirname(os.getcwd()) and not os.getcwd()
    reedSolomonProjectDir = os.path.dirname(os.getcwd()) 
sys.path.insert(0, reedSolomonProjectDir)
from arithmetic import polynomial, gf256

# Constant
SYMBOL_LENGTH = 8
SEED = 7134066
PARITY_LENGTH = 15

import numpy as np
localPRNG = np.random.RandomState(SEED)
CODEWORD_LENGTH = (2 ** SYMBOL_LENGTH) - 1



# Create a generator polynomial for encoding messages, capable of decoding up to designedDistance = 15 erasures
alpha = gf256([0,0,0,0,0,0,1,0])
beta = gf256([0,0,0,0,0,0,1,0])
wan = gf256([0,0,0,0,0,0,0,1])
zro = gf256([0,0,0,0,0,0,0,0])
# Initialize the generator polynomial as (X - alpha)
generatorPolynomialX = polynomial([wan, alpha])
for i in range(PARITY_LENGTH):
    beta = beta * alpha
    linearFactor = polynomial([wan, beta])
    generatorPolynomialX = generatorPolynomialX * linearFactor

# Print the coefficients of the generator polynomial
generatorPolynomialX.printValues()

# generate some data
data = localPRNG.randint(0, 2, (CODEWORD_LENGTH - PARITY_LENGTH) * SYMBOL_LENGTH)
# parse the binary data into 8-bit-long-symbols and cast into (a list of) gf256 symbols
dataParsed = [gf256(data[i * 8 : (i + 1) * 8]) for i in range(CODEWORD_LENGTH  - PARITY_LENGTH)]
# Create a codeword polynomial that contains the data as the coefficients of the higher-than-PARITY_LENGTH powers
codewordX = polynomial(dataParsed).lift(PARITY_LENGTH)
codewordX.printValues()
parityX = codewordX.modulu(generatorPolynomialX)
codewordX.coefficients[0 : PARITY_LENGTH] = parityX.coefficients

# Test - codewordX modulu generator polynomial should now be 0:
generatorPolynomialX.printValues()
codewordX.printValues()
modulusX = codewordX.modulu(generatorPolynomialX)
modulusX.printValues()



