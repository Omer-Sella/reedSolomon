# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 11:41:12 2024

@author: Omer
"""
from keyEquationSolver import *
from arithmetic import binaryFieldElement as galoisElement
from arithmetic import gf128, generateExponentAndLogTables, polynomial
from arithmetic import polynomial as polynomialClass
import numpy as np
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


    #r = np.array([0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0]) #This is equivalent to r(x) = x + x^3 + x^8
    syndromes = np.array([ 1, 1, 1 , 0, 1 , 0, 0])
    #syndromes = np.array([ 0,0,1,0,1,1,1])
    syndromesAsGFElements = []
    for s in syndromes:
        syndromesAsGFElements.append(galoisElement(s))
    cX = keyEquationSolver(polynomialClass, galoisElement, syndromesAsGFElements)
    cX.printValues()
    #assert(cX.coefficients[0].value == 1)
    #assert(cX.coefficients[1].value == 1)
    #assert(cX.coefficients[2].value == 0)
    #assert(cX.coefficients[3].value == 1)
    
    assert(cX.getCoefficient(0).value == 1)
    assert(cX.getCoefficient(1).value == 1)
    assert(cX.getCoefficient(2).value == 0)
    assert(cX.getCoefficient(3).value == 1)


def test_keyEquationSolver_bug_connection_polynomial_for_one_error_order_check():   
    
    #Testing the Todd K. Moon version using example on page 281 (laboratory 6 ex. 2)
    # As well as table 6.5 on page 259
    #|K | S[k]                       | discrepancyK    | cX                         | L | pX   | l | oldDiscrepancy == d_m in Moon
    #|1 | [0,0,0,0,0,0,1]            | [0,0,0,0,0,0,1] | 1 + x                      | 1 | 1    | 1 | 1
    #|2 | \alpha = [0,0,0,0,0,1,0]   | [0,0,0,0,0,1,1] | 1 + [0,0,0,0,0,1,1]x       | 1 | 1    | 2 | 1
    #|3 | \alpha^2 = [0,0,0,0,1,0,0] | 0               | 1 + [0,0,0,0,0,1,0]x       | 1 | 1    | 3 | 1
    
    
    
    #|4 | \alpha^3 = [0,0,0,1,0,0,0] | 0               | 1 + [0,0,0,0,0,1,0]x       | 1 | 1    | 4 | 1
    #|5 | \alpha^4 = [0,0,1,0,0,0,0] | 0               | 1 + [0,0,0,0,0,1,0]x       | 1 | 1    | 5 | 1
    
    syndromes = np.array([[0,0,0,0,0,0,1] , [0,0,0,0,0,1,0], [0,0,0,0,1,0,0],  [0,0,0,1,0,0,0],  [0,0,1,0,0,0,0]])
    syndromeAsGf128 = list(map(gf128, syndromes))
    cX = keyEquationSolver(polynomialClass, gf128, syndromeAsGf128)
    #cX.printValues()
    assert cX.order() == 1
    
def test_keyEquationSolver_bug_connection_polynomial_for_one_error_explicit_calculation():   
    eD, _ =  generateExponentAndLogTables()
    syndromeAsGf128 = []
    for i in range(16):
        syndromeAsGf128.append(gf128(eD[i]))
    cX = keyEquationSolver(polynomialClass, gf128, syndromeAsGf128)
    lambda0 = gf128(1)
    lambda1 = gf128([0,0,0,0,0,1,0])
    #lambda2 = (syndromeAsGf128[3] + (syndromeAsGf128[1] * syndromeAsGf128[1] * syndromeAsGf128[1]) / syndromeAsGf128[0])
    explicitConnectionX = polynomial(coefficients = [lambda1, lambda0])
    assert explicitConnectionX == cX
    
