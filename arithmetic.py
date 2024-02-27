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
    
    def isZero(self):
        return (self.value == 0)

class polynomial():
    # From draft dj_d0p1:
    #— Define c as a binary vector of length 126, and c(x) a polynomial of degree 125 with the coefficients 
    #defined by c (where the bit 0 of c represents the coefficient of power 125).
    #Meaning coefficients[0] is the highest power
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
    
    def ordrer(self):
        #find first non zero:
        order = len(self.coefficients) - 1
        i = 0
        while self.coefficients[i].isZero():
            i = i + 1
            order = order - 1
        return ( order )
    
    def truncate(self):
        while self.order != (len(self.coefficients) - 1):
            self.coefficients.pop(0)
        return self
            
    def plus(self, other):
        newCoefficients = []
        # if other.order() > self.order():
        #     small = self
        #     big = other
        # else:
        #     small = other
        #     big = other
        if other.order() > self.order():
            for i in range(other.order() - self.order()):
                newCoefficients.append(binaryFieldElement(0))
            for i in range(self.order()):
                newCoefficients.append(other.coefficients[i + self.order()])
        if self.order() > other.order():
            for i in range(self.order() - other.order()):
                newCoefficients.append(binaryFieldElement(0))
            for i in range(self.order() - other.order()):
                newCoefficients.append(other.coefficients[i + other.order()])
        return polynomial(coefficients = newCoefficients)
    

    def lift(self, liftBy):
        for i in range(liftBy):
            self.coefficients.append(binaryFieldElement(0))
        return self
    
    def times(self, other):
        newCoefficients = []
        for i in range(len(other.coefficients) , 0 ,-1):
            fieldElement = other.coefficients[i]
            temp = polynomial(self.coefficients)
            temp.lift(len(other.coefficients) - i)
            for j in range(len(temp.coefficients)):
                temp.coefficients[j] = temp.coefficients[j].times(fieldElement)
            return polynomial(coefficients = newCoefficients)        
                
def test_constructor():
    p1 = polynomial([1,0,0])
    assert(p1.order == 3)
    p2 = polynomial([0,1,0,0])
    assert(p2.order == 3)
    return 'OK'

def test_lift():
    p1 = polynomial([1,1,1])
    p1 = p1.lift(10)
    assert (p1.order == (10 + 3 - 1))
    
def test_truncate():
    p0 = polynomial([1,0,0])
    p1 = polynomial([0,1,0,0])
    p1.truncate()
    assert np.all(p0.coefficients == p1.coefficients)
    return 'OK'