# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:36:27 2024

@author: Omer
"""
import numpy as np
import copy
class binaryFieldElement:
    def __init__(self, value):
        if value == 0 or value == 1:
            self.value = value
        else:
            raise("A binary field value can be either 0 or 1 in this class.")
    
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
    def __init__(self, coefficients = None):
        newCoefficients = []
        if coefficients is not None:
            if type(coefficients[0]) == type(binaryFieldElement(value = 0)):
                for i in range(len(coefficients)):
                    newCoefficients.append(copy.deepcopy(coefficients[i]))
            elif np.isscalar(coefficients[0]):
                for i in range(len(coefficients)):
                    newCoefficients.append(binaryFieldElement(coefficients[i]))
            else:
                raise('Çoefficients must be finite field elements or scalalrs')    
        self.coefficients = newCoefficients
    
    def order(self):
        #find first non zero:
        order = len(self.coefficients) - 1
        i = 0
        while self.coefficients[i].isZero():
            i = i + 1
            order = order - 1
        return order 
    
    def truncate(self):
        while ((self.coefficients[0].value == 0) and len(self.coefficients) > 1) : #(len(self.coefficients) - 1) > self.order():
            self.coefficients.pop(0)
        return self
            
    def plus(self, other):
        #newCoefficients = copy.deepcopy(self.coefficients)
        if other.order() > self.order():
            small = self
            big = other
            bigOrder = other.order()
            smallOrder = self.order()
            
        else:
            #print("Other one is smaller")
            small = other
            big = self
            bigOrder = self.order()
            smallOrder = other.order()
        #print("Small order is " + str(smallOrder))
        #print("Big order is " + str(bigOrder))
        newCoefficients = copy.deepcopy(big.coefficients)
        #print(len(newCoefficients))
        for i in range(smallOrder + 1):
            print("now in index " + str(i))
            a = newCoefficients[len(newCoefficients) - 1 - i]
            b = small.coefficients[len(small.coefficients)- 1 - i]
            newCoefficients[bigOrder - i] = a.plus(b)
            print(newCoefficients[bigOrder - i].value)

        
        return polynomial(coefficients = newCoefficients)
    
    def minus(self, other):
        return self.plus(other)
    

    def lift(self, liftBy):
        for i in range(liftBy):
            self.coefficients.append(binaryFieldElement(0))
        return self
    
    def timesScalar(self, gfScalar):
        newCoefficients = copy.deepcopy(self.coefficients)
        for j in range(len(self.coefficients)):
            newCoefficients[j] = newCoefficients[j].times(gfScalar)
        return polynomial(newCoefficients)
    
    def times(self, other):
        newCoefficients = []
        for i in range(len(other.coefficients) , 0 ,-1):
            fieldElement = other.coefficients[i]
            temp = polynomial(self.coefficients)
            temp.lift(len(other.coefficients) - i)
            temp.timesScalar(fieldElement)        
        return polynomial(coefficients = newCoefficients)        
                
    def __eq__(self, other):
        # This is a super lazy implementation, since I actually deepcopy and test equality on the truncated polynomials
        pSelf = copy.deepcopy(self).truncate()
        pOther = copy.deepcopy(other).truncate()
        isEqual = True
        if pSelf.order() == pOther.order():
            for i in range(len(pSelf.coefficients)):
                if pSelf.coefficients[i] != pOther.coefficients[i]:
                    isEqual = False
        else:
            isEqual = False
        
        return isEqual
    
    def printValues(self):
        for e in self.coefficients:
            print(e.value)