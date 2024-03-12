# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:03:01 2024

@author: Omer
"""
from arithmetic import binaryFieldElement as galoisElement
from arithmetic import polynomial as polynomialClass

def keyEquationSolver(polynomialClass, galoisElement, syndromes):
     # Key Equation Solver over an extension binary field
     # See Todd K. Moon page 259
     L = 0 #Current Length of the LFSR
     cX = polynomialClass(coefficients = [1]) #Connection polynomial
     pX = polynomialClass(coefficients = [1]) #The connection polynomial before last change
     l = 1
     oldDiscrepancy = galoisElement(value = 1) #was 1  but I think might be wrong
     for k in range(len(syndromes)):
         discrepancy = galoisElement(value = 0)
         for i in range(L):
             discrepancy = discrepancy + cX[i].times((syndromes[k-i]))
         discrepancy = syndromes + discrepancy
         if discrepancy.value() == 0:
             l = l + 1
         else:
             if ((2 * L) >= k):
                 cX = cX.minus(pX.timesScalar((discrepancy.times(oldDiscrepancy.inverse()))))
                 l = l + 1
             else:
                 tX = cX
                 cX = cX.minus(pX.timesScalar((discrepancy.times(oldDiscrepancy.inverse()))))
                 L = k - L
                 pX = tX
                 oldDiscrepancy = discrepancy
                 l = 1
     return cX

        
def berlekampMassey():
    #An implementation following Richard E. Blahut page 187
    L = 0
    r = 0
    lambdaX = 1
    bX = 1
    return lambdaX

def test_keyEquationSolver():
    import numpy as np
    #r = np.array([0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0]) #This is equivalent to r(x) = x + x^3 + x^8
    syndromes = np.array([ 1, 1, 1 , 0, 1 , 0, 0])
    cX = keyEquationSolver(polynomialClass, galoisElement, syndromes)
    assert(np.all(cX.coefficients == np.array([1,1,0,1])))