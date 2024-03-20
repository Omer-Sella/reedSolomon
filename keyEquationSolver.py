# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:03:01 2024

@author: Omer
"""
from arithmetic import binaryFieldElement as galoisElement
from arithmetic import polynomial as polynomialClass

def keyEquationSolver(polynomialClass, galoisElement, syndromes):
     # Key Equation Solver over an extension binary field
     # See Todd K. Moon page 258
     L = 0 #Current Length of the LFSR
     cX = polynomialClass(coefficients = [1]) #Connection polynomial
     pX = polynomialClass(coefficients = [1]) #The connection polynomial before last change
     l = 1
     oldDiscrepancy = galoisElement(value = 1) #was 1  but I think might be wrong
     for k in range(len(syndromes)):
         K = k + 1
         discrepancy = galoisElement(value = syndromes[k])
         for i in range(L):
             print("i == " + str(i))
             discrepancy = discrepancy + cX.coefficients[i].times((galoisElement(syndromes[K-i])))
         # ???discrepancy = syndromes + discrepancy
         if discrepancy.isZero():
             l = l + 1
         else:
             if (2 * L) >= K: #Note the fix k replaced with K
                 cX = cX - (pX.timesScalar((discrepancy.times(oldDiscrepancy.inverse())))).lift(l)
                 l = l + 1
             else:
                 tX = cX
                 cX = cX - (pX.timesScalar((discrepancy.times(oldDiscrepancy.inverse())))).lift(l)
                 L = K - L
                 pX = tX
                 oldDiscrepancy = discrepancy
                 l = 1
         print("| K | syndromes[k] | disc | l  | L   | old discrepancy |")
         print("| %d | %d            | %d    | %d  | %d   | %d              |" % (K, syndromes[k], discrepancy.value, l, L,  oldDiscrepancy.value))
         
         print("c(x) == ")
         cX.printValues()
         print("p(x) == ")
         pX.printValues()
         
        
     return cX

        
def berlekampMassey():
    #An implementation following Richard E. Blahut page 187
    L = 0
    r = 0
    lambdaX = 1
    bX = 1
    return lambdaX

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
    cX = keyEquationSolver(polynomialClass, galoisElement, syndromes)
    cX.printValues()
    assert(np.all(cX.coefficients == np.array([1,1,0,1])))
