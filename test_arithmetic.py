# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 12:05:03 2024

@author: Omer
"""
from arithmetic import *

def test_fieldPlus():
    a = binaryFieldElement(0)
    b = binaryFieldElement(1)
    c = a.plus(b)
    assert c.value == 1
    c = b.plus(a)
    assert c.value == 1
    

def test_constructor():
    p1 = polynomial([1,0,0])
    assert(p1.order() == 2)
    p2 = polynomial([0,1,0,0])
    assert(p2.order() == 2)
    

def test_lift():
    p1 = polynomial([1,1,1])
    p1 = p1.lift(10)
    assert (p1.order() == (10 + 3 - 1))
    
    
def test_truncate():
    p0 = polynomial([1,0,0])
    p1 = polynomial([0,1,0,0])
    p1 = p1.truncate()
    assert (p0.coefficients[0] == p1.coefficients[0])
    assert (p0.coefficients[1] == p1.coefficients[1])
    assert (p0.coefficients[2] == p1.coefficients[2])
    

def test_plus():
    p0 = polynomial([1,0,0])
    p1 = polynomial([0,1,1])
    p4 = p0.plus(p1)
    p2 = polynomial([1,1,1])
    assert(p0.plus(p1).coefficients[0] == 1)
    assert(np.all(p0.plus(p1).coefficients[1] == 1))
    assert(np.all(p0.plus(p1).coefficients[2] == 1))
    assert(np.all(p1.plus(p2).coefficients[0] == 1))
    assert(np.all(p1.plus(p2).coefficients[1] == 2))
    assert(np.all(p1.plus(p2).coefficients[2] == 2))
    assert(np.all(p2.plus(p0).coefficients[0] == 2))
    assert(np.all(p2.plus(p0).coefficients[1] == 1))
    assert(np.all(p2.plus(p0).coefficients[2] == 1))
    

def test_timesScalar():
    
    coefficients = np.random.randint(0,2,100)
    coefficientsAsBinary = []
    zro = binaryFieldElement(0)
    for e in coefficients:
        coefficientsAsBinary.append(binaryFieldElement(value = e))
    p0 = polynomial(coefficientsAsBinary)
    scalar0 = binaryFieldElement(0)
    scalar1 = binaryFieldElement(1)
    assert (np.all(p0.timesScalar(scalar1).coefficients == p0.coefficients))
    assert (np.all(p0.timesScalar(scalar0).coefficients == zro))
    

def test_mul():
    #x^7 + x^3 + 1
    a = polynomial([1,0,0,0,1,0,0,1])
    #x^9 + x^7 + x^5 + x^4 + x + 1
    b = polynomial([1,0,1,0,1,1,0,0,1,1])
    c = (a * b)
    c.coefficients = c.coefficients %2
    assert(c.order() == 16)
    assert(c.coefficients[16-0] == 1)
    assert(c.coefficients[16-1] == 1)
    assert(c.coefficients[16-2] == 0)
    assert(c.coefficients[16-3] == 1)
    assert(c.coefficients[16-4] == 0)
    assert(c.coefficients[16-5] == 1)
    assert(c.coefficients[16-6] == 0)
    assert(c.coefficients[16-7] == 1)
    assert(c.coefficients[16-8] == 0)
    assert(c.coefficients[16-9] == 1)
    assert(c.coefficients[16-10] == 1)
    assert(c.coefficients[16-11] == 1)
    assert(c.coefficients[16-12] == 0)
    assert(c.coefficients[16-13] == 0)
    assert(c.coefficients[16-14] == 1)
    assert(c.coefficients[16-15] == 0)
    assert(c.coefficients[16-16] == 1)
    

def test_modulu():
    #x^7 + x^3 + 1
    a = polynomial([1,0,0,0,1,0,0,1])
    #x^9 + x^7 + x^5 + x^4 + x + 1
    b = polynomial([1,0,1,0,1,1,0,0,1,1])
    c = polynomial([1,0,1,0,0,1,1,1,0,1,0,1,0,1,0,1,1])
    d = c.modulu(b)
    e = polynomial([0])
    assert (d == e)
    

def test_typeConsistencyGalois128():
    a = gf128(1)
    b = gf128(0)
    c = a.mul(b)
    assert c.__class__ == a.__class__
    

def test_plus_bug1():
    #0*X^16 0*X^15 0*X^14 0*X^13 1*X^12 0*X^11 1*X^10 1*X^9 1*X^8 0*X^7 0*X^6 1*X^5 0*X^4 1*X^3 0*X^2 1*X^1 1*X^0
    a = polynomial(coefficients = [0, 0, 0, 0, 1, 0, 1, 1, 1 ,0, 0 ,1, 0, 1, 0, 1, 1])
    #1*X^12 0*X^11 1*X^10 0*X^9 1*X^8 1*X^7 0*X^6 0*X^5 1*X^4 1*X^3 0*X^2 0*X^1 0*X^0
    b = polynomial(coefficients = [1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0])
    c = a + b 
    c.coefficients = c.coefficients %2
    #The bug was: c = 0*X^16 1*X^15 1*X^14 1*X^13 1*X^12 1*X^11 0*X^10 1*X^9 0*X^8 0*X^7 0*X^6 1*X^5 1*X^4 1*X^3 0*X^2 1*X^1 1*X^0
    testResult = polynomial([1,0,1,0,1,1,0,0,1,1])
    assert c == testResult
    

def test_plus_bug2():
    #0*X^0
    a = polynomial(coefficients=[0])
    #0*X^12 0*X^11 0*X^10 0*X^9 0*X^8 0*X^7 0*X^6 0*X^5 0*X^4 0*X^3 0*X^2 0*X^1 0*X^0
    b = polynomial(coefficients=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    c = a + b
    c.coefficients = c.coefficients %2
    testResult = polynomial([0])
    assert c == testResult
    

def test_evaluate_at_value_bug():
    a = gf128(1)
    p = polynomial(coefficients=[a])
    eD, _ =  generateExponentAndLogTables()
    p.at(gf128(eD[0]))
    
def test_noZeroDivision():
    from itertools import combinations
    eD, _ =  generateExponentAndLogTables()
    for combination in combinations(range(len(eD.keys())), 2):
        assert ((gf128(eD[combination[0]]) * gf128(eD[combination[1]])) != 0)
        
def test_plusCommutativity():
    from itertools import combinations
    eD, _ =  generateExponentAndLogTables()
    for combination in combinations(range(len(eD.keys())), 2):
        a = gf128(eD[combination[0]]) + gf128(eD[combination[1]])
        b = gf128(eD[combination[1]]) + gf128(eD[combination[0]])
        assert (a == b)
        
def test_timesCommutativity():
    from itertools import combinations
    eD, _ =  generateExponentAndLogTables()
    for combination in combinations(range(len(eD.keys())), 2):
        a = gf128(eD[combination[0]]) * gf128(eD[combination[1]])
        b = gf128(eD[combination[1]]) * gf128(eD[combination[0]])
        assert (a == b)
    
# def test_polynomialMultiplication():
#     eD, _ =  generateExponentAndLogTables()
#     generatorX = polynomial(coefficients = [gf128(1)])
#     for i in range(16):
#         linearFactor = polynomial(coefficients = [gf128(1), gf128(eD[i+1])])
#         generatorX = generatorX * linearFactor
#         generatorX.printValues()

def test_chienSearch():
    eD, _ =  generateExponentAndLogTables()
    # g(x) =      x16 +  x14 +    x11 + x10 + x9 + x7 + x5 +  x3 +  x+ 1
    gX = polynomial(coefficients = list(map(gf128,np.array([1,0,1,0,0,1,1,1,0,1,0,1,0,1,0,1,1]))))
    p1X = polynomial(coefficients = list(map(gf128,np.array([1,0,0,0,1,0,0,1]))))
    p2X = polynomial(coefficients = list(map(gf128,np.array([1,0,0,0,1,1,1,1]))))
    #[[63, [0, 0, 0, 1, 0, 1, 0]],
    # [95, [0, 1, 1, 1, 1, 0, 0]],
    # [111, [1, 1, 1, 1, 1, 1, 0]],
    # [119, [1, 1, 0, 0, 0, 1, 0]],
    # [123, [1, 0, 0, 1, 1, 0, 0]],
    # [125, [0, 1, 0, 0, 0, 1, 0]],
    # [126, [1, 0, 0, 0, 1, 0, 0]]]
    
    #[[31, [0, 0, 0, 0, 0, 1, 1]],
    #[62, [0, 0, 0, 0, 1, 0, 1]],
    #[79, [0, 0, 1, 0, 1, 0, 1]],    
    #[103, [0, 0, 0, 0, 1, 1, 1]],
    #[115, [0, 0, 1, 0, 1, 1, 1]],
    #[121, [0, 0, 1, 0, 0, 1, 1]],
    #[124, [0, 0, 1, 0, 0, 0, 1]]]
    p3X = polynomial(coefficients = list(map(gf128,np.array([1,1]))))
    assert (gX == p1X * p2X * p3X * p3X)
    
def test_logExponentRoundtrip():
    a = gf128(1)
    for i in a.exponentTable.keys():
        c = gf128(a.exponentTable[i])
        assert (c.getLog() == i)


def test_exponentiation():
    a = gf128(1)
    b = gf128([0,0,0,0,0,1,0])
    c = gf128([0,0,0,0,0,1,0])
    for i in range(2, 140, 1):
        c = c * b    
        d = gf128(a.exponentTable[(i * b.getLog()) % len(a.exponentTable)])
        assert d == c
    
def test_polynomialRealNumbers():
    a = polynomial(coefficients = [1,0,-1]) # x^2 - 1 
    c = a.at(1)
    assert c == 0
    c = a.at(-1)
    assert c == 0
    
def test_polynomialCatchAlignmentBug():
    a = polynomial(coefficients = [1,0,0]) # x^2
    assert(a.at(2) == 4)
    assert(a.at(3) == 9)
    
def test_polynomialCatchAlignmentBugOverGF128():
    a = polynomial(coefficients = [gf128([0,0,0,0,0,0,1]), gf128([0,0,0,0,0,1,0])]) #x + \alpha
    assert(a.at(gf128([0,0,0,0,0,1,0])) == gf128(0))