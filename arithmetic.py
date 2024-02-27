# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:36:27 2024

@author: Omer
"""
import numpy as np
class binaryFieldElement:
    def __init__(self, value = None):
        if value == 0 or value == 1:
            self.value = value
        else:
            raise Exception("A binary field value can be either 0 or 1 in this class.")
    
    def times(self, other):
        return binaryFieldElement(self.value * other.value)
    
    def inverse(self):
        if self.value == 0:
            raise("Division bu zero is not permitted")
        return binaryFieldElement(value = self.value)
    
    def plus(self, other):
        return binaryFieldElement( ((self.value + other.value) % 2))
    
    def minus(self, other):
        return self.plus(other)

class polynomial():
    def __init__(self, coefficients = None, order = None):
        if coefficients is None:
            newCoefficients = []
        else:
            for i in range(len(coefficients)):
                newCoefficients.append(binaryFieldElement(coefficients[i]))
        if order is not None:
            if order <= len(coefficients):
                raise Exception("Order of a polynomial has to be greater or equal number of coefficients.")
            else:
                for i in range(len(coefficients)):
                    newCoefficients.append(binaryFieldElement(coefficients[i]))
                for i in range(order + 1 - len(coefficients)):
                    newCoefficients.append(binaryFieldElement(0))
        self.coefficients = newCoefficients
    
    def ordre(self):
        return (len(self.coefficients) + 1)
    
    def plus(self, other):
        biggerOrder = self.order()
        newCoefficients = []
        if other.order() > self.order():
            newCoefficients = copy.deepcopy(other.coefficients)
            for i in range(other.order() - self.order()):
                newCoefficients.append(other.coefficients[i + self.order()])
        if self.order() > other.order():
            newCoefficients = copy.deepcopy(self.coefficients)
            for i in range(self.order() - other.order()):
                newCoefficients.append(other.coefficients[i + other.order()])
        return polynomial(coefficients = newCoefficients)
    
    def times(self, other):
        newCoefficients = []
        
    return polynomial(coefficients = newCoefficients)
            
        
            
    def times(self, other):
                