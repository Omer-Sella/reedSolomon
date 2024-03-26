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
    return 'OK'

def test_lift():
    p1 = polynomial([1,1,1])
    p1 = p1.lift(10)
    assert (p1.order() == (10 + 3 - 1))
    
def test_truncate():
    p0 = polynomial([1,0,0])
    p1 = polynomial([0,1,0,0])
    p1 = p1.truncate()
    assert (p0.coefficients[0].value == p1.coefficients[0].value)
    assert (p0.coefficients[1].value == p1.coefficients[1].value)
    assert (p0.coefficients[2].value == p1.coefficients[2].value)
    return 'OK'

def test_plus():
    p0 = polynomial([1,0,0])
    p1 = polynomial([0,1,1])
    p4 = p0.plus(p1)
    p2 = polynomial([1,1,1])
    assert(p0.plus(p1).coefficients[0].value == 1)
    assert(np.all(p0.plus(p1).coefficients[1].value == 1))
    assert(np.all(p0.plus(p1).coefficients[2].value == 1))
    assert(np.all(p1.plus(p2).coefficients[0].value == 1))
    assert(np.all(p1.plus(p2).coefficients[1].value == 0))
    assert(np.all(p1.plus(p2).coefficients[2].value == 0))
    assert(np.all(p2.plus(p0).coefficients[0].value == 0))
    assert(np.all(p2.plus(p0).coefficients[1].value == 1))
    assert(np.all(p2.plus(p0).coefficients[2].value == 1))
    return 'OK'

def test_timesScalar():
    p0 = polynomial(np.random.randint(0,2,100))
    scalar1 = binaryFieldElement(0)
    scalar2 = binaryFieldElement(1)
    assert (np.all(p0.timesScalar(scalar1).coefficients == p0.coefficients))
    assert (np.all(p0.timesScalar(scalar2).coefficients == (0 * p0.coefficients)))
    return 'OK'

def test_mul():
    #x^7 + x^3 + 1
    a = polynomial([1,0,0,0,1,0,0,1])
    #x^9 + x^7 + x^5 + x^4 + x + 1
    b = polynomial([1,0,1,0,1,1,0,0,1,1])
    c = a * b
    assert(c.order() == 16)
    assert(c.coefficients[16-0].value == 1)
    assert(c.coefficients[16-1].value == 1)
    assert(c.coefficients[16-2].value == 0)
    assert(c.coefficients[16-3].value == 1)
    assert(c.coefficients[16-4].value == 0)
    assert(c.coefficients[16-5].value == 1)
    assert(c.coefficients[16-6].value == 0)
    assert(c.coefficients[16-7].value == 1)
    assert(c.coefficients[16-8].value == 0)
    assert(c.coefficients[16-9].value == 1)
    assert(c.coefficients[16-10].value == 1)
    assert(c.coefficients[16-11].value == 1)
    assert(c.coefficients[16-12].value == 0)
    assert(c.coefficients[16-13].value == 0)
    assert(c.coefficients[16-14].value == 1)
    assert(c.coefficients[16-15].value == 0)
    assert(c.coefficients[16-16].value == 1)
    return 'OK'

def test_modulu():
    #x^7 + x^3 + 1
    a = polynomial([1,0,0,0,1,0,0,1])
    #x^9 + x^7 + x^5 + x^4 + x + 1
    b = polynomial([1,0,1,0,1,1,0,0,1,1])
    c = polynomial([1,0,1,0,0,1,1,1,0,1,0,1,0,1,0,1,1])
    d = c.modulu(b)
    e = polynomial([0])
    assert (d == e)
    return 'OK'