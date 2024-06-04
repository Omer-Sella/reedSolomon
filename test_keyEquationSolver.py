# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 11:41:12 2024

@author: Omer
"""
from keyEquationSolver import *
from arithmetic import binaryFieldElement as galoisElement
from arithmetic import polynomial as polynomialClass

def test_keyEquationSolver():
    #Testing the Todd K. Moon version using example on page 281 (laboratory 6 ex. 2)
    # As well as table 6.5 on page 259
    #|K | S[k] | discrepancyK | cX          | L | pX   | l | oldDiscrepancy
    #|1 | 1    | 1            | 1+x         | 1 | 1    | 1 | 1
    #|2 | 1    | 0            | 1+x         | 1 | 1    | 1 | 1
    #|3 | 1    | 0            | 1+x         | 1 | 1    | 1 | 1
    #|4 | 0    | 1            | 1+x+x^3     | 3 | 1+x  | 1 | 1
    #|5 | 1    | 0            | 1+x+x^3     | 3 | 1+x  | 1 | 1
    #|6 | 0    | 0            | 1+x+x^3     | 3 | 1+x  | 1 | 1
    #|7 | 0    | 0            | 1+x+x^3     | 3 | 1+x  | 1 | 1

    import numpy as np
    #r = np.array([0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0]) #This is equivalent to r(x) = x + x^3 + x^8
    syndromes = np.array([ 1, 1, 1 , 0, 1 , 0, 0])
    #syndromes = np.array([ 0,0,1,0,1,1,1])
    syndromesAsGFElements = []
    for s in syndromes:
        syndromesAsGFElements.append(galoisElement(s))
    cX = keyEquationSolver(polynomialClass, galoisElement, syndromesAsGFElements)
    #cX.printValues()
    assert(cX.coefficients[0].value == 1)
    assert(cX.coefficients[1].value == 1)
    assert(cX.coefficients[2].value == 0)
    assert(cX.coefficients[3].value == 1)
    return 'OK'
