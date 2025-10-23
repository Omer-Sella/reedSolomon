# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:36:27 2024

@author: Omer
"""
from logging import setLogRecordFactory
import os, sys
reedSolomonProjectDir = os.environ.get('REEDSOLOMON')
if reedSolomonProjectDir == None: 
     reedSolomonProjectDir = os.getcwd()#"c:/users/omer/reedSolomon/reedSolomon/"
sys.path.insert(0, reedSolomonProjectDir)
import numpy as np
import copy
from abc import ABC, abstractmethod


"""
This file contains some ad-hoc arithmetic over the binary field.
"""
IEEE_BINARY_DTYPE = np.int32


G_8_0 = [1,0,0,0,1,1,1,0,1] #x^8 + x^4 + x^3 + x^2 + 1
G_8_1 = [1,0,0,1,0,1,0,1,1]#x^8 + x^5 + x^3 + x^1 + 1
#x^8 + x^6 + x^4 + x^3 + x^2 + x^1 + 1
#x^8 + x^6 + x^5 + x^1 + 1
#x^8 + x^6 + x^5 + x^2 + 1
#x^8 + x^6 + x^5 + x^3 + 1
#x^8 + x^7 + x^6 + x^1 + 1
#x^8 + x^7 + x^6 + x^5 + x^2 + x^1 + 1

G_7_0 =  [1,0,0,0,0,0,1,1]   #x^7 + x^1 + 1 
G_7_1 =  [1,0,0,0,1,0,0,1] #x^7 + x^3 + 1 
G_7_2 =  [1,0,0,0,1,1,1,1] #x^7 + x^3 + x^2 + x^1 + 1 
G_7_3 =  [1,0,0,1,1,1,0,1] #x^7 + x^4 + x^3 + x^2 + 1 
G_7_4 =  [1,0,1,1,1,1,1,1] #x^7 + x^5 + x^4 + x^3 + x^2 + x^1 + 1 
G_7_5 =  [1,1,0,0,1,0,1,1] #  #x^7 + x^6 + x^3 + x^1 + 1
G_7_6 =  [1,1,0,1,0,1,0,1] #x^7 + x^6 + x^4 + x^2 + 1
G_7_7 =  [1,1,1,0,0,1,0,1] # x^7 + x^6 + x^5 + x^2 + 1
G_7_8 =  [1,1,1,1,0,1,1,1] #x^7 + x^6 + x^5 + x^4 + x^2 + x^1 + 1 


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
    
    def mul(self, other):
        return self.times(other)
    
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
        #         raise('Coefficients must be finite field elements or scalalrs')    
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
    
    def symbolicPlus(self, other):
        if len(other.coefficients) > len(self.coefficients):
            newCoefficients = copy.deepcopy(other.coefficients)
            diff = len(other.coefficients) - len(self.coefficients)
            for i in range(len(self.coefficients)):
                newCoefficients[diff + i] = "((" + newCoefficients[diff + i] + " + " + self.coefficients[i] +") %2)"
        elif len(other.coefficients) < len(self.coefficients):
            newCoefficients = copy.deepcopy(self.coefficients)
            diff = len(self.coefficients) - len(other.coefficients)
            for i in range(len(other.coefficients)):
                newCoefficients[diff + i] = "((" + newCoefficients[diff + i] + " + " + other.coefficients[i] + ")%2)"
        else:
            newCoefficients = copy.deepcopy(self.coefficients)
            diff = 0
            for i in range(len(other.coefficients)):
                newCoefficients[diff + i] = "((" + newCoefficients[diff + i] + " + " + other.coefficients[i] + ")%2)"
    
    def minus(self, other):
        return self.plus(other)
    
    def symbolicMinus(self, other):
        return self.symbolicPlus(other)
    

    def lift(self, liftBy):
        if liftBy < 0:
            raise ValueError("Lifting a polynomial can only be done for positive integers")
        if liftBy > 0:
            newCoefficients = list(self.coefficients)
            for i in range(liftBy):
                newCoefficients.append(self.coefficients[0].__class__(0)) 
            newCoefficients = np.array(newCoefficients)    
            self.coefficients = newCoefficients
        return self
    
    def timesScalar(self, gfScalar):
        """
        Note that this function does not alter self, instead it returns a new polynomial which is the result.
        """
        newCoefficients = copy.deepcopy(self.coefficients)
        for j in range(len(self.coefficients)):
            newCoefficients[j] = newCoefficients[j] * gfScalar
        return polynomial(newCoefficients)
    
    def symbolicTimesScalar(self, string):
        newCoefficients = copy.deepcopy(self.coefficients)
        for j in range(len(self.coefficients)):
            newCoefficients[j] = "((" + newCoefficients[j] + " * " + string + ") %2)"
        return polynomial(newCoefficients)
    
    def times(self, other):
        i = 0
        length = len(other.coefficients)
        result = polynomial(coefficients = [self.coefficients[0].__class__(0)])#[0]) # BUG !!! this should be polynomial(self.coefficients[0].__class__(0)), but the change triggers a bug in exponent table generation
        while i < length: 
            fieldElement = other.coefficients[i]
            temp = polynomial(self.coefficients)
            temp.lift(length - 1 - i)
            temp = temp.timesScalar(fieldElement)    
            #print(f"result is {result.coefficients} and temp is {temp.coefficients}")
            result = result.plus(temp)
            #print(f"After adding them the result is: {result.coefficients}")
            i = i + 1
        return result

    def symbolicTimes(self, other):
        i = 0
        length = len(other.coefficients)
        result = polynomial(coefficients = ['0'])#self.coefficients[0].__class__(0)])
        while i < length: 
            fieldElement = other.coefficients[i]
            temp = polynomial(self.coefficients)
            temp.lift(length - 1 - i)
            temp = temp.symbolicTimesScalar(fieldElement)    
            result = result.symbolicPlus(temp)
            i = i + 1
        return result

    def polynomialMultiplication(self, other):
        i = 0
        length = len(other.coefficients)
        result = polynomial(coefficients = [self.coefficients[0].__class__(0)]) # BUG !!! this should be polynomial(self.coefficients[0].__class__(0)), but the change triggers a bug in exponent table generation
        while i < length: 
            fieldElement = other.coefficients[i]
            temp = polynomial(self.coefficients)
            temp.lift(length - 1 - i)
            temp = temp.timesScalar(fieldElement)    
            result = result.plus(temp)
            i = i + 1
        return result

    def ignoreFromDegree(self, degreeToIgnoreFrom):
        #if (degreeToTruncateFrom + 1) >= len(self.coefficients):
        if (degreeToIgnoreFrom ) <= self.order(): 
            #newCoefficients = self.coefficients[0 : (degreeToTruncateFrom + 1)]
            # Recall that the leading coefficient is the FIRST non zero coefficient, and if the order of a polynomial is D, then there are at least D+1 coefficients (with possibly leading 0s)
            newCoefficients = self.coefficients[len(self.coefficients) - self.getLeadingCoefficientIndex() - degreeToIgnoreFrom : ]
            return polynomial(coefficients = newCoefficients)
        else:
            return self
    
    def scalarType(self):
        return self.coefficeints[0].__class__
        
    def d(self):
        # This only works over GF(2), since (d/dx)(x^2) == 0 over GF(2)
        if self.order() == 0:
            return polynomial(self.scalarType()(0))
        else:
            newCoefficients = []
            # reduce the degree of every monomial by 1
            for m in range(1, len(self.coefficients), 1):
                print(m)
                if m % 2 == 0:
                    # Indices with even power will have zero contribution to the derivative over GF(2)
                    newCoefficients.append(self.scalarType(0))
                else:
                    newCoefficients.append(self.coefficients[m])
            result = polynomial(newCoefficients)
            return result
    
        
    def modulu(self, divisor):
        # Safety - no division by zero
        if divisor.isZero():
            raise ValueError("Divisor cannot be the 0 polynomial.")
        # Safety - the leading coefficient of the divisor is not zero
        if divisor.coefficients[0] == 0:
            raise ValueError("An attempt was made to divide by a polynomial which leading coefficient is 0. Truncate the divisor before dividing.")
        one = self.coefficients[0].__class__(1)
        remainder = polynomial(copy.deepcopy(self.coefficients))
        divisorOrder = divisor.order()
        remainderOrder = remainder.order()
        divisorLeadingCoefficientInverse = one / divisor.coefficients[0]
        #divisorLeadingCoefficientInverse = one // divisor.coefficients[0]
        if divisorOrder == 0:
            remainder = polynomial( coefficients = [0])
        while remainderOrder >= divisorOrder and not remainder.isZero():
            i = remainder.getLeadingCoefficientIndex()
            leadingCoefficient = remainder.coefficients[i]#self.coefficients[i]
            #kill ther leading coefficient
            leadingCoefficientKiller = leadingCoefficient * divisorLeadingCoefficientInverse
            temp = polynomial(copy.deepcopy(divisor.coefficients))
            temp.lift(remainderOrder - divisorOrder)
            # Important !!! Note that while pX.timesScalar(scalar) returns the desired result, it does not alter pX !!!
            temp = temp.timesScalar(leadingCoefficientKiller)
            remainder = remainder + temp #temp.timesScalar(leadingCoefficientKiller) # Note + rather than -, all over a binary field
            if np.isscalar(self.coefficients[0]):
                remainder.coefficients = remainder.coefficients %2
            remainderOrder = remainder.order()
        # Omer: I don't think the line below is necessary anymore, since remainder is already declared as a polynomial.    
        remainder = polynomial(remainder.coefficients[-len(divisor.coefficients) + 1 : ])
        return remainder

    def symbolicModulu(self, divisor):
        # WARNING - the assumption here is that the coefficients are only 0 or 1 !!!
        # This function is really only used when calculating symbolic multiplication of gf(2**M) elements
        remainder = polynomial(copy.deepcopy(self.coefficients))
        divisorOrder = divisor.order()
        symbolicOrder = remainder.order()
        divisorSymbolicOrder = divisor.order()
        i = 0 #remainder.getLeadingCoefficientIndex()
        while symbolicOrder >= divisorOrder and not remainder.isZero():
            
            #kill ther leading coefficient
            fieldElementInverse = remainder.coefficients[i] #.inverse()
            temp = polynomial(copy.deepcopy(divisor.coefficients))
            temp.lift(symbolicOrder - divisorSymbolicOrder)
            temp.symbolicTimesScalar(fieldElementInverse)
            remainder = remainder.symbolicPlus(temp) 
            #if np.isscalar(self.coefficients[0]):
            #    remainder.coefficients = remainder.coefficients %2
            symbolicOrder = symbolicOrder - 1
            remainder.coefficients.popleft()
        remainder = polynomial(remainder.coefficients[-len(divisor.coefficients) + 1 : ])
        return remainder
    
    def at(self, evaluationPoint):
        # No safety ! The multiplication between the evaluation point and the coefficients needs to make sense.
        # Initialize result as the zero of galois field  of the same class as evaluationPoint 
        result = evaluationPoint.__class__(0)
        if hasattr(evaluationPoint, 'logTable'):
            if evaluationPoint == 0:
                result = copy.deepcopy(self.coefficients[-1])
            else:
                logEvaluationPoint = evaluationPoint.getLog()
                exponentArray = (np.arange(len(self.coefficients), -1 , -1) * logEvaluationPoint) % len(evaluationPoint.exponentTable)
                # I could have vectorized the addition and multiplication but let's see if the log and exponent alone are enoughj
                #elementWiseMultiply = np.array([evaluationPoint.__class__(evaluationPoint.exponentTable[i]) * self.coefficients[i] for i in range(len(self.coefficients))])
                for i in range(len(self.coefficients)):
                    result = result + (evaluationPoint.__class__(evaluationPoint.exponentTable[exponentArray[i]]) * self.coefficients[i])
                
        else: 
            # Use the brute force method
            
            # Initialize the helper gfElement to be the 1 of galois field  of the same class as evaluationPoint 
            powerOfEvaluationPoint = evaluationPoint.__class__(1)
            # Polynomials are leading coefficient at index 0
            for i in range(len(self.coefficients)):    
                result = result + (self.getCoefficient(i) * powerOfEvaluationPoint)
                powerOfEvaluationPoint = powerOfEvaluationPoint * evaluationPoint
        return result
    
    def getCoefficient(self, index):
        return self.coefficients[len(self.coefficients) - (index + 1)]
    
    def asArray(self):   
        if np.isscalar(self.coefficients[0]):
            return self.coefficients
        else:
            return np.asarray([c.asArray() for c in self.coefficients]).flatten()
    
    
                
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
                return np.all(pSelf.coefficients == pOther.coefficients)
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
            

class gfBase(polynomial, ABC):
    """
    A common API for all manifestations of galois fields over the binary field.
    The user must specify:
    polynomialName : one of the constants representing polynomials at the top of this file.
    pathToInverseTable : we could actually do without this, but I'm keeping it for now for speedup
    pathToLogTable : we could actually do without this, but I'm keeping it for now for speedup
    pathToExponentTable : we could actually do without this, but I'm keeping it for now for speedup
    
    """
    
    @property
    @abstractmethod
    def polynomialName(self):
        raise NotImplementedError
    
    @property
    @abstractmethod
    def lengthInBits(self):
        # Needs to be len(polynomialName)
        raise NotImplementedError    
    
    @property
    @abstractmethod
    def pathToInverseTable(self):
        raise NotImplementedError
    
    @property
    @abstractmethod
    def pathToLogTable(self):
        raise NotImplementedError
    
    @property
    @abstractmethod
    def pathToExponentTable(self):
        raise NotImplementedError
    
    @property
    @abstractmethod
    def pathToTimesTable(self):
        pass
    
    @property
    @classmethod
    def inverseTable(cls):
        return np.load(cls.pathToInverseTable, allow_pickle = True).item()  
        
    
    @property
    @classmethod
    def exponentTable(cls):
        return np.load(cls.pathToExponentTable, allow_pickle = True).item()
    
    @property
    @classmethod
    def logTable(cls):
        return np.load(cls.pathToLogTable, allow_pickle = True).item()
    
    @property
    @classmethod
    def timesTable(cls):
        if cls.pathToTimesTable is not None:
            return np.load(cls.pathToTimesTable, allow_pickle = True).item()
        else:
            return None
    
    # It is assumed that the polynomial used in the child class is one of the named constant polynomials at the top of this file
    @property
    @classmethod
    def loadGeneratorPolynomial(cls):
        cls.generatorPolynomial = polynomial(cls.polynomialName)
    
    
    def __init__(self, value):
        
        if hasattr(value, '__len__'):
            if len(value) == self.lengthInBits:
                super().__init__(coefficients = value)
            else:
                #print("Class of provided value is " + str(value.__class__))
                raise(f"An element in GF({2 ** self.lengthInBits}) is a {self.lengthInBits}-tuple of binary values. Please avoid ambiguity by stating all {self.lengthInBits} coefficients. ")
        elif np.isscalar(value) and (value == 0 or value == 1):
            coefficients = np.zeros(self.lengthInBits, IEEE_BINARY_DTYPE)
            coefficients[-1] = value
            # WARNING ! make sure polynomial is the first class in the inheritance list
            super().__init__(coefficients = coefficients)
        else:
            raise(f"Class of provided value is {value.__class__}. An element in GF({2 ** self.lengthInBits} is a {self.lengthInBits}-tuple of binary values. Please avoid ambiguity by stating all {self.lengthInBits} coefficients.")
    
    def mul(self, other):        
        tempResult = self.times(other)
        tempResult = tempResult.modulu(self.generatorPolynomial)
        result = self.__class__(value = tempResult.coefficients)
        return result
    
    def inverse(self):
        return self.__class__(self.inverseTable[str(self.getValue())])
    
    def __mul__(self, other):
        return self.__class__(self.timesTable["".join(str(e) for e in self.getValue())]["".join(str(c) for c in other.getValue())])
    
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
        return self * other.inverse()

    def __add__(self, other):
        return self.plus(other)
     
    def __sub__(self, other):
        return self.minus(other)
      
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
    
    def getLog(self):
        if self == 0 :
            raise ValueError('Log is not defined for 0')
        else:
            key = ''.join(map(str, self.coefficients))
        return self.logTable[key]
    
class gf256Generic(gfBase):
    @property
    def polynomialName(self):
        return G_8_0
    @property
    def lengthInBits(self):
        return 8
        
    @property
    def pathToInverseTable(self):
        return reedSolomonProjectDir + "/cachedArithmetic/gf256Inverse_G_8_0.npy"
        
    @property
    def pathToExponentTable(self):
        return reedSolomonProjectDir + "/cachedArithmetic/gf256Exponent_G_8_0.npy"
        
    @property
    def pathToLogTable(self):
        return reedSolomonProjectDir + "/cachedArithmetic/gf256Log_G_8_0.npy"
        
    @property
    def pathToTimesTable(self):
        return reedSolomonProjectDir + "/cachedArithmetic/gf256TimesTable_G_8_0.npy"
        
        
        
class gf128Generic(gfBase):
    super
    def polynomialName(self):
        return G_7_1
    @property
    def lengthInBits(self):
        return 7
        
    @property
    def pathToInverseTable(self):
        return reedSolomonProjectDir + "/cachedArithmetic/gf128Inverse.npy"
        
    @property
    def pathToExponentTable(self):
        return reedSolomonProjectDir + "/cachedArithmetic/gf128Exponent.npy"
        
    @property
    def pathToLogTable(self):
        return reedSolomonProjectDir + "/cachedArithmetic/gf128Log.npy"
        
    @property
    def pathToTimesTable(self):
        return None


class gf128(polynomial):
    """
    Ad-hoc implementation of gf128 arithmetic, since this arithmetic class has several possible specific optimizations 
    """
    
    pathToInverseTable = reedSolomonProjectDir + "/cachedArithmetic/gf128Inverse.npy"
    pathToExponentTable = reedSolomonProjectDir + "/cachedArithmetic/gf128Exponent.npy"
    pathToLogTable = reedSolomonProjectDir + "/cachedArithmetic/gf128Log.npy"
    inverseTable = np.load(pathToInverseTable, allow_pickle = True).item()
    exponentTable = np.load(pathToExponentTable, allow_pickle = True).item()
    logTable = np.load(pathToLogTable, allow_pickle = True).item()
    generatorPolynomial = polynomial([1,0,0,0,1,0,0,1])

    
    
    def __init__(self, value):
        if hasattr(value, '__len__'):
            if len(value) == 7:
                super().__init__(coefficients = value)
            else:
                print("Class of provided value is " + str(value.__class__))
                raise("An element in GF(128) is a 7-tuple of binary values. Please avoid ambiguity by stating all 7 coefficients. ")
        elif np.isscalar(value) and (value == 0 or value == 1):
            coefficients = np.zeros(7, IEEE_BINARY_DTYPE)
            coefficients[-1] = value
            super().__init__(coefficients = coefficients)
        else:
            raise("Class of provided value is " + str(value.__class__) + "An element in GF(128) is a 7-tuple of binary values. Please avoid ambiguity by stating all 7 coefficients.")
    
    def mul(self, other):
        ########### GF128 mul
        tempResult = self.times(other)
        tempResult.coefficients = tempResult.coefficients %2
        # Omer Sella: note that there was a choice of polynomila here, namely: p(x) = x^7 + x^3 + 1
        tempResult = tempResult.modulu(self.generatorPolynomial)#polynomial([1,0,0,0,1,0,0,1]))
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
    
    def getLog(self):
        if self == 0 :
            raise ValueError('Log is not defined for 0')
        else:
            key = ''.join(map(str, self.coefficients))
        return self.logTable[key]


class gf256(polynomial):
    
    polynomialName = "G_8_0"
    pathToInverseTable = reedSolomonProjectDir + "/cachedArithmetic/gf256Inverse_" + polynomialName + ".npy"
    pathToLogTable = reedSolomonProjectDir + "/cachedArithmetic/gf256Log_" + polynomialName + ".npy"
    pathToExponentTable = reedSolomonProjectDir + "/cachedArithmetic/gf256Exponent_"+ polynomialName + ".npy"
    inverseTable = np.load(pathToInverseTable, allow_pickle = True).item()  
    exponentTable = np.load(pathToExponentTable, allow_pickle = True).item()
    logTable = np.load(pathToLogTable, allow_pickle = True).item()
    pathToTimesTable = reedSolomonProjectDir + "/cachedArithmetic/gf256TimesTable_" + polynomialName + ".npy"
    timesTable = np.load(pathToTimesTable, allow_pickle = True).item()
    # Irreducible polynomials from https://www.partow.net/programming/polynomials/index.html#deg08
    #In 177-1 it seems like they used the polynomial x^8 + x^7 + 0 + x^5 + x^4 + 0 + 0 + x + 1 which is not irreducible / primitive 
    generatorPolynomial = polynomial(G_8_0)
    def __init__(self, value):
        if hasattr(value, '__len__'):
            if len(value) == 8:
                super().__init__(coefficients = value)
            else:
                print("Class of provided value is " + str(value.__class__))
                raise("An element in GF(256) is a 8-tuple of binary values. Please avoid ambiguity by stating all 8 coefficients. ")
        elif np.isscalar(value) and (value == 0 or value == 1):
            coefficients = np.zeros(8, IEEE_BINARY_DTYPE)
            coefficients[-1] = value
            super().__init__(coefficients = coefficients)
        else:
            raise("Class of provided value is " + str(value.__class__) + "An element in GF(256) is a 8-tuple of binary values. Please avoid ambiguity by stating all 8 coefficients.")
    
    def mul(self, other):
        ########### GF256 mul
        tempResult = self.times(other)
        tempResult.coefficients = tempResult.coefficients %2
        tempResult = tempResult.modulu(self.generatorPolynomial)
        result = self.__class__(value = tempResult.coefficients)
        return result
        
    def inverse(self):
        if np.all(self.coefficients == 0):
            raise ZeroDivisionError
        return self.__class__(self.inverseTable[str(self.getValue())]) 

    
    def __mul__(self, other):
        answer = self.__class__(self.timesTable["".join(str(e) for e in self.getValue())]["".join(str(c) for c in other.getValue())])
        return answer
    
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
        return self * other.inverse()

    def __add__(self, other):
        return self.plus(other)
     
    def __sub__(self, other):
        return self.minus(other)
     
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
    
    def getLog(self):
        if self == 0 :
            raise ValueError('Log is not defined for 0')
        else:
            key = ''.join(map(str, self.coefficients))
        return self.logTable[key]


def generateTimesTable(gfType = gf128, generatorAsList = [0,0,0,0,0,1,0]):
    a = gfType(generatorAsList)
    b = gfType(generatorAsList)
    generator = gfType(generatorAsList)
    timesTable = dict()
    zro = gfType([0] * len(generatorAsList))
    # Populate zero times b and b times zero
    timesTable["".join(str(c) for c in zro.getValue())] = dict()
    timesTable["".join(str(c) for c in zro.getValue())]["".join(str(e) for e in zro.getValue())] = zro.getValue()
    for i in range(2 ** len(generatorAsList) - 1):        
        timesTable["".join(str(c) for c in zro.getValue())]["".join(str(e) for e in b.getValue())] = zro.getValue()
        timesTable["".join(str(c) for c in b.getValue())] = dict()
        timesTable["".join(str(c) for c in b.getValue())]["".join(str(e) for e in zro.getValue())] = zro.getValue()
        b = b.mul(generator)
    # Reset b
    b = gfType(generatorAsList)
    for i in range(2 ** len(generatorAsList) - 1):
        for j in range(2 ** len(generatorAsList) - 1):
            result = a.mul(b)
            timesTable["".join(str(e) for e in a.getValue())]["".join(str(ee) for ee in b.getValue())]  = result.getValue()
            #for coefficient in result.coefficients:
                #assert(coefficient == 0 or coefficient == 1)
            b = b.mul(generator)
        a = a.mul(generator)
    return timesTable

def generateExponentAndLogTables(gfType = gf128, generatorAsList = [0,0,0,0,0,1,0]):
    exponentTable={}
    logarithmTable={}
    a = gfType(generatorAsList) #gf128([0,0,0,0,0,1,0])
    wanAsList = [0] * len(generatorAsList)
    wanAsList[len(wanAsList) - 1] = 1
    b = gfType(wanAsList) #gf128([0,0,0,0,0,0,1])
    wanAsString = "".join([str(e) for e in wanAsList])
    f = []
    stringF = '' 
    for e in a.coefficients:
        f.append(e)
        stringF = stringF + str(e)
    exponentTable[0] = wanAsList #[0,0,0,0,0,0,1]
    #exponentTable[1] = f
    logarithmTable[wanAsString] = 0
    logarithmTable[stringF] = 1
    for i in range(1, (2 ** len(generatorAsList)) - 1, 1):
        b = b * a
        f = []
        stringF = ''
        for e in b.coefficients:
            f.append(e)
            stringF = stringF + str(e)
        exponentTable[i] = f
        logarithmTable[stringF] = i
    return exponentTable, logarithmTable

def generateInverseTable(gfType = gf128, generatorAsList = [0,0,0,0,0,1,0]):
    # A pretty lazy implementation of finding an inverse, we're only using this once in a lifetime, so simple and readable.
    inverseDictionary = {}
    wanAsList = [0] * len(generatorAsList)
    wanAsList[len(wanAsList) - 1] = 1
    a = gfType(wanAsList) #gf128([0,0,0,0,0,0,1])
    temp = gfType(wanAsList) #gf128([0,0,0,0,0,0,1])
    key = str(a.getValue())
    b = gfType(generatorAsList) #gf128([0,0,0,0,0,1,0])
    ONE = gfType(wanAsList) #gf128([0,0,0,0,0,0,1])
    inverseDictionary[key] = temp.getValue()
    for i in range(1,(2 ** len(generatorAsList) ) - 1, 1):
        print(i)
        a = a * b
        key = str(a.getValue())
        temp = gfType(wanAsList)# gf128([0,0,0,0,0,1,0])
        while temp * a != ONE:
            temp = temp * b
        inverseDictionary[key] = temp.getValue()
    return inverseDictionary