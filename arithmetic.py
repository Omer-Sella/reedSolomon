# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:36:27 2024

@author: Omer
"""
import numpy as np
import copy
"""
This file contains some ad-hoc arithmetic over the binary field.
"""
IEEE_BINARY_DTYPE = np.int32

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
    
    def __add__(self, other):
        return self.plus(other)
    
    def __sub__(self, other):
        return self.minus(other)
    
    def getValue(self):
        return self.value

class polynomial():
    # From draft dj_d0p1:
    #— Define c as a binary vector of length 126, and c(x) a polynomial of degree 125 with the coefficients 
    #defined by c (where the bit 0 of c represents the coefficient of power 125).
    #Meaning coefficients[0] is the highest power
    def __init__(self, coefficients = None):
        newCoefficients = []
        if np.isscalar(coefficients[0]):
            for i in range(len(coefficients)):
                newCoefficients.append(binaryFieldElement(coefficients[i]))
        else:
            newCoefficients = copy.deepcopy(coefficients)
        self.coefficients = newCoefficients
        # if coefficients is not None:
        #     if type(coefficients[0]) == type(binaryFieldElement(value = 0)):
        #         for i in range(len(coefficients)):
        #             newCoefficients.append(copy.deepcopy(coefficients[i]))
        #     elif np.isscalar(coefficients[0]):
        #         for i in range(len(coefficients)):
        #             newCoefficients.append(binaryFieldElement(coefficients[i]))
        #     else:
        #         raise('Çoefficients must be finite field elements or scalalrs')    
        # self.coefficients = newCoefficients

    def isZero(self):
        zro = True
        for coefficient in self.coefficients:
            if not coefficient.isZero():
                zro = False
        return zro
        
    def order(self):
        #find first non zero:
        length = len(self.coefficients)
        order = length - 1
        i = 0
        while (i < length) and self.coefficients[i].isZero():
            i = i + 1
            order = order - 1
        return order
    
    def getLeadingCoefficientIndex(self):
        i = 0
        length = len(self.coefficients)
        while ( i < length) and self.coefficients[i].isZero():
            i = i + 1
        return i
        
    
    def truncate(self):
        while ((self.coefficients[0].value == 0) and len(self.coefficients) > 1) : #(len(self.coefficients) - 1) > self.order():
            self.coefficients.pop(0)
        return self
            
    def plus(self, other):
        if other.order() > self.order():
            small = self
            big = other
            bigOrder = other.order()
            smallOrder = self.order()
            
        else:
            small = other
            big = self
            bigOrder = self.order()
            smallOrder = other.order()
        newCoefficients = copy.deepcopy(big.coefficients)
        for i in range(smallOrder + 1):
            a = newCoefficients[len(newCoefficients) - 1 - i]
            b = small.coefficients[len(small.coefficients)- 1 - i]
            newCoefficients[bigOrder - i] = a.plus(b)
        return polynomial(coefficients = newCoefficients)
    
    def minus(self, other):
        return self.plus(other)
    

    def lift(self, liftBy):
        assert(liftBy >=0 )
        for i in range(liftBy):
            self.coefficients.append(binaryFieldElement(0))
        return
    
    def timesScalar(self, gfScalar):
        newCoefficients = []#copy.deepcopy(self.coefficients)
        for j in range(len(self.coefficients)):
            newCoefficients.append(self.coefficients[j].times(gfScalar))
        return polynomial(newCoefficients)
    
    def times(self, other):
        i = 0
        length = len(other.coefficients)
        result = polynomial([0])
        while i < length: 
            fieldElement = other.coefficients[i]
            temp = polynomial(self.coefficients)
            temp.lift(length - 1 - i)
            temp = temp.timesScalar(fieldElement)    
            result = result.plus(temp)
            i = i + 1
        return result

    def ignoreFromDegree(self, degreeToTruncateFrom):
        newCoefficients = self.coefficients[0 : (degreeToTruncateFrom + 1)]
        return polynomial(coefficients = newCoefficients)
        
    def d(self):
        # This only works over GF(2):
        if len(self.coefficients) % 2 == 0:
            length = len(self.coefficients)
        else:
            length = len(self.coefficients) + 1
        mask = np.arange(0, length, 2)
        newCoefficients = copy.deepcopy(self.coefficients)
        # Indices with even power will have zero contribution to the derivative over GF(2)
        newCoefficients[mask] = 0
        #Now reduce the degree of every mono by 1
        newCoefficients.pop()
        return polynomial(newCoefficients)
    
        
    def modulu(self, divisor):
        assert (not divisor.isZero())
        remainder = polynomial(copy.deepcopy(self.coefficients))
        while remainder.order() >= divisor.order() and not remainder.isZero():
            i = remainder.getLeadingCoefficientIndex()
            #if not remainder.coefficients[0].isZero():
            #kill ther leading coefficient
            fieldElementInverse = remainder.coefficients[i].inverse()
            temp = polynomial(copy.deepcopy(divisor.coefficients))
            temp.lift(remainder.order() - divisor.order())
            temp.timesScalar(fieldElementInverse)
            remainder = remainder + temp
        remainder = polynomial(remainder.coefficients[-len(divisor.coefficients) + 1 : ])
        return remainder
    
    def at(self, evaluationPoint):
        # Initialize result as the zero of galois field  of the same class as evaluationPoint 
        result = evaluationPoint.__class__(0)
        # Initialize the helper gfElement to be the 1 of galois field  of the same class as evaluationPoint 
        gfElement = evaluationPoint.__class__(1)
        for i in range(len(self.coefficients)):
            temp = self.coefficients[i].getValue()
            print(temp)
            helper = gfElement.mul(evaluationPoint.__class__(temp))
            print(helper.__class__)
            result = result + ( helper )
            gfElement = gfElement.mul(evaluationPoint)
        return result
    
    
                
    def __eq__(self, other):
        # This is a super lazy implementation, since I actually deepcopy and test equality on the truncated polynomials
        pSelf = copy.deepcopy(self).truncate()
        pOther = copy.deepcopy(other).truncate()
        isEqual = True
        if pSelf.order() == pOther.order():
            for i in range(len(pSelf.coefficients)):
                if pSelf.coefficients[i].value != pOther.coefficients[i].value:
                    isEqual = False
        else:
            isEqual = False
        
        return isEqual
    
    def __add__(self, other):
        return(self.plus(other))
    
    def __sub__(self, other):
        # Addition and sbtraction modulu 2 are both bitwise xor
        return(self.plus(other))
    
    def __mul__(self, other):
        return self.times(other)
    
    def printValues(self):
        for element in self.coefficients:
            print(element.value)
            
class gf128(polynomial):
    def __init__(self, value):
        if hasattr(value, '__len__'):
            if len(value) == 7:
                super().__init__(coefficients = value)
            else:
                print("Class of provided value is " + str(value.__class__))
                raise("An element in GF(128) is a 7-tuple of binary values. Please avoid ambiguity by stating all 7 coefficients. ")
        elif value == 0 or value == 1:
            coefficients = np.zeros(7, IEEE_BINARY_DTYPE)
            coefficients[-1] = value
            super().__init__(coefficients = coefficients)
        else:
            print("Class of provided value is " + str(value.__class__))
            raise("An element in GF(128) is a 7-tuple of binary values. Please avoid ambiguity by stating all 7 coefficients.")
    
    def mul(self, other):
        tempResult = self.times(other)
        tempResult = tempResult.modulu(polynomial([1,0,0,0,1,0,0,1]))
        print(tempResult.coefficients)
        result = gf128(value = tempResult.coefficients)
        return result
        
    def inverse(self):
        return
    
    def binaryMul(self, other):
        if other.value == 0:
            return gf128(value = 0)
        else:
            return gf128(value = self.coefficients)
    
    def getValue(self):
        return self.coefficients

def generateExponentAndLogTables():
    exponentTable={}
    logarithmTable={}
    a = gf128([0,0,0,0,0,1,0])
    b = gf128([0,0,0,0,0,1,0])
    f = []
    stringF = '' 
    for e in b.coefficients:
        f.append(e.value)
        stringF = stringF + str(e.value)
    exponentTable[0] = [0,0,0,0,0,0,1]
    exponentTable[1] = f
    logarithmTable[stringF] = 1
    for i in range(2,127,1):
        b = b * a
        b = b.modulu(polynomial([1,0,0,0,1,0,0,1]))
        f = []
        stringF = '' 
        for e in b.coefficients:
            f.append(e.value)
            stringF = stringF + str(e.value)
        exponentTable[i] = f
        logarithmTable[stringF] = i
    return exponentTable, logarithmTable


