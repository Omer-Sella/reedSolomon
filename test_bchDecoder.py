# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:18:16 2024

@author: Omer
"""
import os, sys
projectDir = os.environ.get('IEEE8032DJ')
#You don't have to define an environment variable, but then you have to give a path to the project here
if projectDir == None: 
     projectDir = "c:/Users/Omer/802.3/"
sys.path.insert(0, projectDir)
reedSolomonProjectDir = os.environ.get('REEDSOLOMON')
if True: #reedSolomonProjectDir == None: 
     reedSolomonProjectDir = "c:/users/omer/reedSolomon/reedSolomon/"
sys.path.insert(0, reedSolomonProjectDir)
from arithmetic import binaryFieldElement as galoisElement
from arithmetic import polynomial as polynomialClass
from arithmetic import generateExponentAndLogTables, gf128, generateExponentAndLogTables
from ieee8023dj_d0p1 import bchEncoder
from bchDecoder import bchDecoder
import numpy as np


def test_bchDecoder():
    zeroData = np.zeros(110)
    encodedZeroData = bchEncoder(zeroData)
    eD, _ =  generateExponentAndLogTables()
    #print("Done generating log and exp dictionaries.")
    correctedVector, correctionVector = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthActual = 126, codewordLengthMaximal = 126)
    assert (np.all(correctedVector == 0))
    
