# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:36:27 2024

@author: Omer
"""
import os, sys
reedSolomonProjectDir = os.environ.get('REEDSOLOMON')
if reedSolomonProjectDir == None: 
     reedSolomonProjectDir = "c:/users/omer/reedSolomon/reedSolomon/"
sys.path.insert(0, reedSolomonProjectDir)
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
    
    def __mul__(self, other):
        return self.times(other)
    
    def __eq__(self, other):
        result = False
        if other.__class__ == self.__class__:
            result = (other.getValue() == self.getValue())
        elif other == 0 or other == 1:
            result = (self.getValue() == other)
        else:
            raise
        return result
    
    def __div__(self, other):
        return self.mul(other.inverse())

    def __truediv__(self, other):
        return self.mul(other.inverse())
    
    def getValue(self):
        return self.value

class polynomial():
    # From draft dj_d0p1:
    #— Define c as a binary vector of length 126, and c(x) a polynomial of degree 125 with the coefficients 
    #defined by c (where the bit 0 of c represents the coefficient of power 125).
    #Meaning coefficients[0] is the highest power
    def __init__(self, coefficients = None):
        #Polynomial is a class over GF(2) 
        #Other than that it should be agnostic to the underlying arithmetic type (an element if GF(2^k) for some k), and uses np.array as a container
        newCoefficients = np.array(coefficients)
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
            if not coefficient == 0:
                zro = False
        return zro
        
    def order(self):
        #find first non zero:
        length = len(self.coefficients)
        order = length - 1
        i = 0
        while (i < length) and (self.coefficients[i] == 0):
            i = i + 1
            order = order - 1
        return order
    
    def getLeadingCoefficientIndex(self):
        i = 0
        length = len(self.coefficients)
        while ( i < length - 1) and (self.coefficients[i] == 0):
            i = i + 1
        return i
        
    
    def truncate(self):
        newCoefficients = self.coefficients[ self.getLeadingCoefficientIndex() : ]
        self.coefficients = newCoefficients
        return self
            
    def old_plus(self, other):
        self.printValues()
        other.printValues()
        if other.order() > self.order():
            small = self
            big = other
        else:
            small = other
            big = self
        bigOrder = self.order()
        smallOrder = other.order()
        newCoefficients = copy.deepcopy(big.coefficients)
        for i in range(len(small.coefficients)):
            a = newCoefficients[ - i -1]
            b = small.coefficients[- i -1]
            newCoefficients[- i -1] = a.plus(b)
        return polynomial(coefficients = newCoefficients)
    
    def plus(self, other):
        if len(other.coefficients) > len(self.coefficients):
            newCoefficients = copy.deepcopy(other.coefficients)
            diff = len(other.coefficients) - len(self.coefficients)
            for i in range(len(self.coefficients)):
                newCoefficients[diff + i] = (newCoefficients[diff + i] + self.coefficients[i])
        elif len(other.coefficients) < len(self.coefficients):
            newCoefficients = copy.deepcopy(self.coefficients)
            diff = len(self.coefficients) - len(other.coefficients)
            for i in range(len(other.coefficients)):
                newCoefficients[diff + i] = (newCoefficients[diff + i] + other.coefficients[i])
        else:
            newCoefficients = copy.deepcopy(self.coefficients)
            diff = 0
            for i in range(len(other.coefficients)):
                newCoefficients[diff + i] = (newCoefficients[diff + i] + other.coefficients[i])
            
        return polynomial(coefficients = newCoefficients)
    
    
    def minus(self, other):
        return self.plus(other)
    

    def lift(self, liftBy):
        assert(liftBy >=0 )
        if liftBy > 0:
            newCoefficients = list(self.coefficients)
            for i in range(liftBy):
                newCoefficients.append(self.coefficients[0].__class__(0)) 
        #newCoefficients = np.zeros((liftBy + len(self.coefficients)), dtype = IEEE_BINARY_DTYPE)
        #newCoefficients[0 : len(self.coefficients)] = self.coefficients
            newCoefficients = np.array(newCoefficients)    
            self.coefficients = newCoefficients
        return self
    
    def timesScalar(self, gfScalar):
        newCoefficients = copy.deepcopy(self.coefficients)
        for j in range(len(self.coefficients)):
            newCoefficients[j] = newCoefficients[j] * gfScalar
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
        # Safety - no division by zero
        assert (not divisor.isZero())
        # Warning - this is only over a binary field !
        one = self.coefficients[0].__class__(1)
        remainder = polynomial(copy.deepcopy(self.coefficients))
        divisorOrder = divisor.order()
        if divisorOrder == 0:
            remainder = polynomial( coefficients = [0])
        while remainder.order() >= divisorOrder and not remainder.isZero():
            i = remainder.getLeadingCoefficientIndex()
            #kill ther leading coefficient
            fieldElementInverse = one / remainder.coefficients[i]
            temp = polynomial(copy.deepcopy(divisor.coefficients))
            temp.lift(remainder.order() - divisor.order())
            temp.timesScalar(fieldElementInverse)
            remainder = remainder + temp # Note + rather than -, all over a binary field
            if np.isscalar(self.coefficients[0]):
                remainder.coefficients = remainder.coefficients %2
            #print("And the new remainder is:")
            #remainder.printValues()
        remainder = polynomial(remainder.coefficients[-len(divisor.coefficients) + 1 : ])
        return remainder
    
    def at(self, evaluationPoint):
        # No safety ! The multiplication between the evaluation point and the coefficients needs to make sense.
        # Initialize result as the zero of galois field  of the same class as evaluationPoint 
        
        result = evaluationPoint.__class__(0)
        # Initialize the helper gfElement to be the 1 of galois field  of the same class as evaluationPoint 
        powerOfEvaluationPoint = evaluationPoint.__class__(1)
        for i in range(len(self.coefficients)):
            temp = self.coefficients[i]
            #temp = temp * gfElement # Note the change - instead of casting temp into the same type as the evaluation point, we assume they are of the same type !
            #temp = (evaluationPoint.__class__(temp)).mul(gfElement)
            temp = temp * powerOfEvaluationPoint
            result = result + temp# * powerOfEvaluationPoint #.plus(temp)
            powerOfEvaluationPoint = powerOfEvaluationPoint * evaluationPoint
        return result
    
    
                
    def __eq__(self, other):
        # This is a super lazy implementation, since I actually deepcopy and test equality on the truncated polynomials
        pSelf = copy.deepcopy(self)
        pSelf = pSelf.truncate()
        pOther = copy.deepcopy(other)
        pOther = pOther.truncate()
        isEqual = True
        if pSelf.order() == pOther.order():
            if np.isscalar(pSelf.coefficients[0]):
                return np.all(pSelf.coefficients == pOther.coefficients)
            else:
                return np.all(pSelf.coefficients.getValue() == pOther.coefficients.getValue())
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
        string = ""
        if np.isscalar(self.coefficients[0]):
            for power in range(len(self.coefficients)):
                string += (" " + str(self.coefficients[power]) + "*X^" + str(len(self.coefficients) - power - 1))
        else:
            for power in range(len(self.coefficients)):
                string += (" " + str(self.coefficients[power].getValue()) + "*X^" + str(len(self.coefficients) - power - 1))
        print(string)
            
class gf128(polynomial):
    
    pathToInverseTable = reedSolomonProjectDir + "/gf128Inverse.npy"
    inverseTable = np.load(pathToInverseTable, allow_pickle = True).item()
    
    def __init__(self, value):
        if hasattr(value, '__len__'):
            if len(value) == 7:
                super().__init__(coefficients = value)
            else:
                print("Class of provided value is " + str(value.__class__))
                raise("An element in GF(128) is a 7-tuple of binary values. Please avoid ambiguity by stating all 7 coefficients. ")
        elif np.isscalar(value) and (value == 0 or value == 1):
            #print("debugging value issue")
            #print(value)
            coefficients = np.zeros(7, IEEE_BINARY_DTYPE)
            coefficients[-1] = value
            super().__init__(coefficients = coefficients)
        else:
            raise("Class of provided value is " + str(value.__class__) + "An element in GF(128) is a 7-tuple of binary values. Please avoid ambiguity by stating all 7 coefficients.")
    
    def mul(self, other):
        #print("Inside mul")
        #print(other.__class__)
        tempResult = self.times(other)
        tempResult = tempResult.modulu(polynomial([1,0,0,0,1,0,0,1]))
        #print(tempResult.coefficients)
        result = gf128(value = tempResult.coefficients)
        return result
        
    def inverse(self):
        if self.inverseTable is not None:
            return gf128(self.inverseTable[str(self.getValue())])
        else:
            # Omer Sella: consider a verilog oriented implementation here, could be slow to run on a CPU but better than not working.
            raise
    
    def binaryMul(self, other):
        if other.value == 0:
            return gf128(value = 0)
        else:
            return gf128(value = self.coefficients)
    
    def getValue(self):
        return self.coefficients
    
    def plus(self, other):
        coefficients = (self.coefficients + other.coefficients) %2
        return gf128(coefficients)
    
    def __div__(self, other):
        return self.mul(other.inverse())
    
    def __truediv__(self, other):
        return self.mul(other.inverse())

    def __add__(self, other):
        return self.plus(other)
     
    def __sub__(self, other):
        return self.minus(other)
     
    def __mul__(self, other):
        return self.mul(other)
     
    def __eq__(self, other):
        result = False
        if other.__class__ == self.__class__:
            result = (np.all(other.getValue() == self.getValue()))
        elif other == 0:
            result = np.all(self.getValue() == 0)
        elif other == 1:
            result = (np.all(self.getValue()[0 : -1] == 0) and (self.getValue()[-1] == 1))
        else:
            raise
        return result
     
    def __ne__(self, other):
        return (not (self == other))

def generateExponentAndLogTables():
    exponentTable={}
    logarithmTable={}
    a = gf128([0,0,0,0,0,1,0])
    b = gf128([0,0,0,0,0,1,0])
    f = []
    stringF = '' 
    for e in b.coefficients:
        f.append(e)
        stringF = stringF + str(e)
    exponentTable[0] = [0,0,0,0,0,0,1]
    exponentTable[1] = f
    logarithmTable[stringF] = 1
    for i in range(2,127,1):
        #print(i)
        b = b * a
        #print(b.coefficients)
        b = b.modulu(polynomial([1,0,0,0,1,0,0,1]))
        f = []
        stringF = '' 
        for e in b.coefficients:
            f.append(e)
            stringF = stringF + str(e)
        exponentTable[i] = f
        logarithmTable[stringF] = i
    return exponentTable, logarithmTable

def generateInverseTable():
    # A pretty lazy implementation of finding an inverse, we're only using this once so why bother.
    inverseDictionary = {}
    a = gf128([0,0,0,0,0,0,1])
    temp = gf128([0,0,0,0,0,0,1])
    key = str(a.getValue())
    b = gf128([0,0,0,0,0,1,0])
    one = gf128([0,0,0,0,0,0,1])
    inverseDictionary[key] = temp.getValue()
    for i in range(1,127,1):
        a = a * b
        key = str(a.getValue())
        temp = gf128([0,0,0,0,0,1,0])
        while temp * a != one:
            temp = temp * b
        inverseDictionary[key] = temp.getValue()
    return inverseDictionary


class gf256(polynomial):
    
    #pathToInverseTable = reedSolomonProjectDir + "/gf256Inverse.npy"
    #inverseTable = np.load(pathToInverseTable, allow_pickle = True).item()
    
    def __init__(self, value):
        if hasattr(value, '__len__'):
            if len(value) == 8:
                super().__init__(coefficients = value)
            else:
                print("Class of provided value is " + str(value.__class__))
                raise("An element in GF(256) is a 8-tuple of binary values. Please avoid ambiguity by stating all 8 coefficients. ")
        elif np.isscalar(value) and (value == 0 or value == 1):
            #print("debugging value issue")
            #print(value)
            coefficients = np.zeros(8, IEEE_BINARY_DTYPE)
            coefficients[-1] = value
            super().__init__(coefficients = coefficients)
        else:
            raise("Class of provided value is " + str(value.__class__) + "An element in GF(256) is a 8-tuple of binary values. Please avoid ambiguity by stating all 8 coefficients.")
    
    def mul(self, other):
        #print("Inside mul")
        #print(other.__class__)
        tempResult = self.times(other)
        tempResult = tempResult.modulu(polynomial([1,0,0,0,1,1,1,0,1])) #x^8 + x^4 + x^3 + x^2 + 1 from Todd K. Moon page 243, example 6.9
        #print(tempResult.coefficients)
        result = self.__class__(value = tempResult.coefficients)
        return result
        
    def inverse(self):
        if self.inverseTable is not None:
            return self.__class__(self.inverseTable[str(self.getValue())])
        else:
            # Omer Sella: consider a verilog oriented implementation here, could be slow to run on a CPU but better than not working.
            raise
    
    def binaryMul(self, other):
        if other.value == 0:
            return self.__class__(value = 0)
        else:
            return self.__class__(value = self.coefficients)
    
    def getValue(self):
        return self.coefficients
    
    def plus(self, other):
        coefficients = (self.coefficients + other.coefficients) %2
        return self.__class__(coefficients)
    
    def __div__(self, other):
        return self.mul(other.inverse())
    
    def __truediv__(self, other):
        return self.mul(other.inverse())

    def __add__(self, other):
        return self.plus(other)
     
    def __sub__(self, other):
        return self.minus(other)
     
    def __mul__(self, other):
        return self.mul(other)
     
    def __eq__(self, other):
        result = False
        if other.__class__ == self.__class__:
            result = (np.all(other.getValue() == self.getValue()))
        elif other == 0:
            result = np.all(self.getValue() == 0)
        elif other == 1:
            result = (np.all(self.getValue()[0 : -1] == 0) and (self.getValue()[-1] == 1))
        else:
            raise
        return result
     
    def __ne__(self, other):
        return (not (self == other))