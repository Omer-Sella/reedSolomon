# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:03:01 2024

@author: Omer
"""
#from arithmetic import binaryFieldElement as galoisElement
#from arithmetic import polynomial as polynomialClass
import copy
from arithmetic import polynomial

def keyEquationSolver(polynomialClass, galoisElementClass, syndromes):
     # Key Equation Solver over an extension binary field
     # See Todd K. Moon page 258
     # No safety - Syndromes need to be passed as a list of syndromes of class galoisElementClass
     # assert (syndromes[0].__class__ == galoisElementClass)
     L = 0 #Current Length of the LFSR
     #cX = polynomialClass(coefficients = [1]) #Connection polynomial
     #pX = polynomialClass(coefficients = [1]) #The connection polynomial before last change
     cX = polynomial(coefficients = [galoisElementClass(1)]) #Connection polynomial
     pX = polynomial(coefficients = [galoisElementClass(1)]) #The connection polynomial before last change
     hX = polynomial(coefficients = [galoisElementClass(0)])
     tX = polynomial(coefficients = [galoisElementClass(0)])
     l = 1
     oldDiscrepancy = galoisElementClass(value = 1) #was 1  but I think might be wrong
     for k in range(len(syndromes)):
         K = k + 1
         discrepancy = syndromes[k]
         if k > 0 :
             i = 1
             while i <= L: #for i in range(1, L, 1):
                 #diff =  cX.coefficients[(i+1)] * ((syndromes[k-(i+1)]))
                 ### BUG ! our notation of polynomials is leading coefficient is actually at index 0 ! That's why I introduced the function "getCoefficient, which returns coefficients[len(coefficients) - (i + 1)]
                 helper = cX.getCoefficient(i) * syndromes[k - i]
                 discrepancy = discrepancy + (helper)
                 i = i + 1
         if discrepancy == 0:
             l = l + 1
         else:
             if (2 * L) >= K: #Note the fix k replaced with K
                 scalar = discrepancy.times(oldDiscrepancy.inverse())
                 print(f"{scalar.getValue()}")
                 hX = pX.timesScalar(scalar)
                 hX.lift(l)
                 print(f"hX is:")
                 hX.printValues()
                 cX = cX - hX
                 l = l + 1
             else:
                 tX = cX
                 # BUG ! discrepancy is gf element, as is oldDiscrepancy.inverse(), but .times() is not a method of gf element (it is a method of polynomial)
                 # Right now the fix is to switch to * instead of .times()
                 #hX = pX.timesScalar(discrepancy.times(oldDiscrepancy.inverse()))
                 scalar = discrepancy * oldDiscrepancy.inverse()
                 print(f"{scalar.getValue()}")
                 hX = pX.timesScalar(scalar)
                 hX.lift(l)
                 print(f"hX is:")
                 hX.printValues()
                 cX = cX - hX
                 L = K - L
                 pX = tX
                 oldDiscrepancy = copy.deepcopy(discrepancy)
                 l = 1
         #print("| K | syndromes[k]    | disc | l  | L   | old discrepancy |")
         #print("| " +str(K) +"| " +str(syndromes[k].getValue()) + "   | "+ str(discrepancy.getValue()) + "  | " + str(l) + "  | " +  str(L) + "   | " + str(oldDiscrepancy.getValue()) +"      |" )
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
    
