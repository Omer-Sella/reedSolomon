# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 12:05:03 2024

@author: Omer
"""
from email import generator
from arithmetic import *
#from symbol import parameters

def test_fieldPlus():
    a = binaryFieldElement(0)
    b = binaryFieldElement(1)
    c = a.plus(b)
    assert c.value == 1
    c = b.plus(a)
    assert c.value == 1
    #return 'OK'
    

def test_constructor():
    p1 = polynomial([1,0,0])
    assert(p1.order() == 2)
    p2 = polynomial([0,1,0,0])
    assert(p2.order() == 2)
    #return 'OK'

def test_lift():
    p1 = polynomial([1,1,1])
    p1 = p1.lift(10)
    assert (p1.order() == (10 + 3 - 1))
    #return 'OK'
    
def test_truncate():
    p0 = polynomial([1,0,0])
    p1 = polynomial([0,1,0,0])
    p1 = p1.truncate()
    assert (p0.coefficients[0] == p1.coefficients[0])
    assert (p0.coefficients[1] == p1.coefficients[1])
    assert (p0.coefficients[2] == p1.coefficients[2])
    #return 'OK'

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
    #return 'OK'

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
    #return 'OK'

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
    #return 'OK'

def test_modulu():
    #x^7 + x^3 + 1
    a = polynomial([1,0,0,0,1,0,0,1])
    #x^9 + x^7 + x^5 + x^4 + x + 1
    b = polynomial([1,0,1,0,1,1,0,0,1,1])
    c = polynomial([1,0,1,0,0,1,1,1,0,1,0,1,0,1,0,1,1])
    d = c.modulu(b)
    e = polynomial([0])
    assert (d == e)
    #return 'OK'

def test_typeConsistencyGalois128():
    a = gf128(1)
    b = gf128(0)
    c = a.mul(b)
    assert c.__class__ == a.__class__
    #return 'OK'

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
    #return 'OK'

def test_plus_bug2():
    #0*X^0
    a = polynomial(coefficients=[0])
    #0*X^12 0*X^11 0*X^10 0*X^9 0*X^8 0*X^7 0*X^6 0*X^5 0*X^4 0*X^3 0*X^2 0*X^1 0*X^0
    b = polynomial(coefficients=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    c = a + b
    c.coefficients = c.coefficients %2
    testResult = polynomial([0])
    assert c == testResult
    #return 'OK'

def test_evaluate_at_value_bug():
    a = gf128(1)
    p = polynomial(coefficients=[a])
    eD, _ =  generateExponentAndLogTables()
    p.at(gf128(eD[0]))
    #return 'OK'

def test_noZeroDivision():
    from itertools import combinations
    #eD, _ =  generateExponentAndLogTables()
    eD = gf128.exponentTable
    for combination in combinations(range(len(eD.keys())), 2):
        assert ((gf128(eD[combination[0]]) * gf128(eD[combination[1]])) != 0)
    eD = gf256.exponentTable
    for combination in combinations(range(len(eD.keys())), 2):
        assert ((gf256(eD[combination[0]]) * gf256(eD[combination[1]])) != 0)
    #return 'OK'

def test_plusCommutativity():
    from itertools import combinations
    #eD, _ =  generateExponentAndLogTables()
    eD = gf128.exponentTable
    for combination in combinations(range(len(eD.keys())), 2):
        a = gf128(eD[combination[0]]) + gf128(eD[combination[1]])
        b = gf128(eD[combination[1]]) + gf128(eD[combination[0]])
        assert (a == b)
    #return 'OK'        
def test_timesCommutativity():
    from itertools import combinations
    eD, _ =  generateExponentAndLogTables()
    for combination in combinations(range(len(eD.keys())), 2):
        a = gf128(eD[combination[0]]) * gf128(eD[combination[1]])
        b = gf128(eD[combination[1]]) * gf128(eD[combination[0]])
        assert (a == b)
    
def test_polynomialMultiplication():
    eD, _ =  generateExponentAndLogTables()
    generatorX = polynomial(coefficients = [gf128(1)])
    for i in range(16):
        linearFactor = polynomial(coefficients = [gf128(1), gf128(eD[i+1])])
        generatorX = generatorX * linearFactor
        generatorX.printValues()

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
        assert (c.getLog() == i )


def test_exponentiation():
    a = gf128(1)
    b = gf128([0,0,0,0,0,1,0])
    c = gf128([0,0,0,0,0,1,0])
    for i in range(2, 140, 1):
        c = c * b    
        d = gf128(a.exponentTable[(i * b.getLog()) % len(a.logTable)])
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
    
def test_reproduceDivisionBug():
    # After investigation, it actually turns out this is a multiplication bug rather than division (multiplication test added).
    aux = gf128([0, 0, 0, 1, 0, 0, 1])
    denominator = gf128([1, 0, 0, 1, 1, 0, 1])
    bug = aux / denominator
    for c in bug.coefficients:
        assert ((c == 0) or (c==1))

def test_reproduceMultiplicationBug():
    aux1 = gf128([0, 0, 0, 1, 0, 0, 1])
    bug = aux1*aux1
    for c in bug.coefficients:
        assert ((c == 0) or (c==1))

def test_gf256CachedMultiplication():
    from itertools import product
    for lhsList in product(range(2), repeat = 8):
        for rhsList in product(range(2), repeat = 8):
            lhs = gf256(lhsList)
            rhs = gf256(rhsList)
            assert lhs * rhs == rhs * lhs
            assert lhs * rhs == lhs.mul(rhs)
    return

def test_moduluOverLargerFields():
    parityLength = 15
    # Create a generator polynomial for encoding messages, capable of decoding up to designedDistance = 15 erasures
    alpha = gf256([0,0,0,0,0,0,1,0])
    beta = gf256([0,0,0,0,0,0,1,0])
    wan = gf256([0,0,0,0,0,0,0,1])
    zro = gf256([0,0,0,0,0,0,0,0])
    # Initialize the generator polynomial as (X - alpha)
    generatorX = polynomial([wan, alpha])
    for i in range(parityLength):
        beta = beta * alpha
        linearFactor = polynomial([wan, beta])
        generatorX = generatorX * linearFactor
    #The degree of generatorX should now be parityLength + 1, and the leading coefficient should be 1    
    assert(generatorX.order() == (parityLength + 1))
    assert(generatorX.coefficients[0] == gf256(1))
    
    messageX = polynomial([gf256(1)]).lift(parityLength + 1)
    parityX = messageX.modulu(generatorX)
    temp = generatorX.ignoreFromDegree(parityLength + 1)
    assert(parityX == temp)

def test_gf256Multiplication():
    wan = gf256(1)
    keys = wan.timesTable.keys()
    for lhs in keys:
        for rhs in keys:
            c = gf256([int(t) for t in lhs])
            d = gf256([int(h) for h in rhs])
            e = c.mul(d)
            for b in e.coefficients:
                assert(b == 0 or b == 1)
    
def test_reproduceGf256MulBug():
    c = gf256([0, 0, 0, 1, 1, 1, 0, 1])
    d = gf256([0, 0, 0, 0, 0, 0, 1, 1])
    e = c.mul(d)
    
    for t in e.coefficients:
        assert(t == 0 or t == 1)

def test_comparisonGf128():
    c = gf128([0, 0, 1, 1, 1, 0, 1])
    d = gf128([0, 0, 0, 0, 0, 1, 1])
    e = c.mul(d)
    
    for t in e.coefficients:
        assert(t == 0 or t == 1)
# def test_reproduceModuluBug():
#     # This was not a bug, but a faulty test !
#     # This test originated from another test where the leading coefficient of the message
#     # was 1, and its order was the same as the generator polynomial, 
#     # in which case the modulus is (exactly) the lower coefficients of the generator polynomial.
    
#     parityLength = 15
#     # Create a generator polynomial for encoding messages, capable of decoding up to designedDistance = 15 erasures
#     alpha = gf256([0,0,0,0,0,0,1,0])
#     beta = gf256([0,0,0,0,0,0,1,0])
#     wan = gf256([0,0,0,0,0,0,0,1])
#     zro = gf256([0,0,0,0,0,0,0,0])
#     # Initialize the generator polynomial as (X - alpha)
#     generatorX = polynomial([wan, alpha])
#     for i in range(parityLength):
#         beta = beta * alpha
#         linearFactor = polynomial([wan, beta])
#         generatorX = generatorX * linearFactor
#     messageX = polynomial([gf256([1,1,1,1,0,1,1,0])]).lift(parityLength + 1)
#     parityX = messageX.modulu(generatorX)
#     temp = generatorX.ignoreFromDegree(parityLength + 1)
#     temp.printValues()
#     parityX.printValues()
#     assert(parityX == temp)
    
def test_polynomialTimesOverLargerFields():
    parityLength = 15
    # Create a generator polynomial for encoding messages, capable of decoding up to designedDistance = 15 erasures
    alpha = gf256([0,0,0,0,0,0,1,0])
    beta = gf256([0,0,0,0,0,0,1,0])
    wan = gf256([0,0,0,0,0,0,0,1])
    zro = gf256([0,0,0,0,0,0,0,0])
    # Initialize the generator polynomial as (X - alpha)
    generatorX = polynomial([wan, alpha])
    for i in range(parityLength):
        beta = beta * alpha
        linearFactor = polynomial([wan, beta])
        generatorX = generatorX * linearFactor
    someElement = gf256([1,1,1,1,0,1,1,0])
    temp = generatorX.timesScalar(someElement)
    #print(temp.coefficients[0].getValue())
    assert(temp.coefficients[0] == someElement)

def test_multiplicationModuluBug():
    # The bug was when gf128(exponentTable[5]) * gf128(exponentTable[6]) was attemted, and the computation never ended when entering the modulu calculation.
    eD = gf128.exponentTable
    assert ((gf128(eD[5]) * gf128(eD[6])) != 0)
    #return 'OK'
    
def test_timesScalarBugReproduction():
    gfScalar = gf256([0, 0, 1, 1, 1, 1, 1, 0])
    coefficients = list(map(gf256, [[0,0, 0, 0, 0, 1, 1, 1], [0, 1, 0, 0, 1, 0, 0, 1], [1, 1, 1, 0, 1,0, 0, 0], [0, 1, 1, 1, 0, 1, 1, 1], [0, 0, 0, 0, 0, 0, 0, 1]]))
    pX = polynomial(coefficients)
    hX = pX.timesScalar(gfScalar)
    
    for c in hX.coefficients:
        for b in c.coefficients:
            assert (b == 0 or b == 1)

def test_gf256MultiplicationTable():
    wan = gf256(1)
    keys = wan.timesTable.keys()
    for lhs in keys:
        for rhs in keys:
            for b in wan.timesTable[lhs][rhs]:
                assert(b == 0 or b == 1)

def test_adHocgf128SameAsUsingBase():
    from arithmetic import reedSolomonProjectDir, gfBase, G_7_1
    class gf128UsingAbstract(gfBase):
        @property
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
    e = gf128UsingAbstract(1)
    
           
if __name__ == "__main__":
    #test_multiplicationModuluBug()
    #test_lift()
    #test_truncate()
    
    #test_noZeroDivision()
    
    #test_chienSearch()
    #test_constructor()
    #test_evaluate_at_value_bug()
    
    #test_logExponentRoundtrip()
    #test_exponentiation()
    #test_modulu()
    #test_polynomialRealNumbers()
    #test_reproduceMultiplicationBug()
    #test_gf256CacheedMultiplication()
    #test_moduluOverLargerFields()
    #test_polynomialTimesOverLargerFields()
    #test_gf256MultiplicationTable()
    #test_gf256Multiplication()
    test_gf256CachedMultiplication()
    #test_comparisonGf128()
    #test_reproduceGf256MulBug()
    #test_timesScalarBugReproduction()
    #test_reproduceModuluBug()
    #test_adHocgf128SameAsUsingBase()
    pass