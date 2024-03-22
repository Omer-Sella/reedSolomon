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
         #print("Disc value == ")
         #print(discrepancy.value)
         for i in range(L):
             #print("k is:" + str(k))
             #print("i is : "+ str(i))
             diff =  cX.coefficients[(i+1)].times((galoisElement(syndromes[k-(i+1)])))
             discrepancy = discrepancy + diff
             #print("Diff value == ")
             #print(diff.value)
         # ???discrepancy = syndromes + discrepancy
         if discrepancy.isZero():
             l = l + 1
         else:
             if (2 * L) >= K: #Note the fix k replaced with K
                 #print("****************")    
                 cX = cX - (pX.timesScalar((discrepancy.times(oldDiscrepancy.inverse())))).lift(l)
                 l = l + 1
             else:
                 
                 tX = cX
                 cX = cX - (pX.timesScalar((discrepancy.times(oldDiscrepancy.inverse())))).lift(l)
                 L = K - L
                 pX = tX
                 oldDiscrepancy = discrepancy
                 l = 1
         #print("| K | syndromes[k] | disc | l  | L   | old discrepancy |")
         #print("| %d | %d            | %d    | %d  | %d   | %d              |" % (K, syndromes[k], discrepancy.value, l, L,  oldDiscrepancy.value))
         
         #print("c(x) == ")
         #cX.printValues()
         #print("p(x) == ")
         #pX.printValues()
         
        
     return cX
 
def findErrorEvaluator(syndromesPolynomial, connectionPolynomial, t):
    omegaX = syndromesPolynomial * connectionPolynomial
    omegaX = omegaX.ignoreFromDegree(2 * t)
    return omegaX


def forney(errorEvaluatorPolynomial, errorConnectionPolynomiol, listOfElements):
    errorConnectionPolynomiolDerivative = errorConnectionPolynomiol.d()
    errorValues = []
    for element in listOfElements:
        # Note that over GF(2) minus and plus are the same, and that the following term should have a minus sign in the non binary case.
        errorValues.append(errorConnectionPolynomiolDerivative.at(element).inverse() * errorEvaluatorPolynomial(element))
    return errorValues
    
