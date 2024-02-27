# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:03:01 2024

@author: Omer
"""
from binaryField import galoisElement, polynomialClass

def keyEquationSolver(polynomialClass, galoisElement, syndromes):
     # Key Equation Solver over an extension binary field
     L = 0
     cX = polynomialClass(coefficients = [galoisElement(value = 1)])
     pX = polynomialClass(coefficients = [galoisElement(value = 1)])
     l = 1
     oldDiscrepancy = galoisElement(value = 1) #was 1  but I think might be wrong
     for k in range(len(syndromes)):
         discrepancy = galoisElement(value = 0)
         for i in range(L):
             discrepancy = discrepancy + cX[i].times(syndromes[k-i])
         discrepancy = syndromes + discrepancy
         if discrepancy.value() == 0:
             l = l + 1
         else:
             if ((2 * L) >= k):
                 cX = cX.minus(pX.times(discrepancy.times(oldDiscrepancy.inverse())))
                 l = l + 1
             else:
                 tX = cX
                 cX = cX.minus(pX.times(discrepancy.times(oldDiscrepancy.inverse())))
                 L = k - L
                 pX = tX
                 oldDiscrepancy = discrepancy
                 l = 1
     return cX

        
    