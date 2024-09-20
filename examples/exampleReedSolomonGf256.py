# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 09:31:34 2024

@author: Omer Sella
"""
import os, sys
reedSolomonProjectDir = os.environ.get('REEDSOLOMON')
if reedSolomonProjectDir == None: 
    #Note that this is an example, so the relevant modules will be in its parent directory, hence os.path.dirname(os.getcwd()) and not os.getcwd()
    reedSolomonProjectDir = os.path.dirname(os.getcwd()) 
sys.path.insert(0, reedSolomonProjectDir)
from arithmetic import polynomial, gf256

### Reed Solomon code over GF(2^8) example

# Create a generator polynomial for encoding messages, capable of decoding up to designedDistance = 15 erasures
designedDistance = 15
alpha = gf256([0,0,0,0,0,0,1,0])
beta = gf256([0,0,0,0,0,0,1,0])
wan = gf256(1)
# Initialize the generator polynomial as (X - alpha)
generatorPolynomial = polynomial([wan, alpha])
for i in range(designedDistance - 1):
    beta = beta * alpha
    linearFactor = polynomial([wan, beta])
    generatorPolynomial = generatorPolynomial * linearFactor

generatorPolynomial.printValues()
