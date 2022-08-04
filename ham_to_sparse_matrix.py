from collections import defaultdict

import scipy
import scipy.sparse

import sympy as sp
import numpy as np

import copy
import math

from timer import Timer

pauliX = np.array([[0,1],[1,0]])
pauliY = np.array([[0,-1j],[1j,0]])
pauliZ = np.array([[1,0],[0,-1]])

def a_element(i,j):
    if(i-1 == j):
        return math.sqrt(i)
    else:
        return 0

def adag_element(i,j):
    if(i+1 == j):
        return math.sqrt(i+1)
    else:
        return 0

def aMat(n):
    return np.array( [[ a_element(i,j) for i in range(n)] for j in range(n)] )

def aDagMat(n):
    return np.array( [[ adag_element(i,j) for i in range(n)] for j in range(n)] )

def xMat():
    return 0.5*(pauliX+complex(0,1)*pauliY)

def xDagMat():
    return 0.5*(pauliX-complex(0,1)*pauliY)



def site_subs(cutoff, Nsites, aops, adags, xs, xdags):
    allSubs = {}
    
    # TODO from 0 to Nsites otherwise error.
    for n in range(-1,Nsites+1): 
        allSubs[aops[n]]=sp.Matrix(aMat(cutoff))
        allSubs[adags[n]]=sp.Matrix(aDagMat(cutoff))
        allSubs[xs[n]]=sp.Matrix(xMat())
        allSubs[xdags[n]]=sp.Matrix(xDagMat())
        
    return allSubs


def site_subs_boson(cutoff, Nsites, aops, adags):
    allSubs = {}
    
    # TODO from 0 to Nsites otherwise error.
    for n in range(-1,Nsites+1): 
        allSubs[aops[n]]=sp.Matrix(aMat(cutoff))
        allSubs[adags[n]]=sp.Matrix(aDagMat(cutoff))
        
    return allSubs



def convert_to_matrix(expr, cutoff, Nsites, aops, adags, xs, xdags):
    # start with a matrix of zeros
    fullHam = scipy.sparse.coo_matrix(np.zeros([(cutoff**Nsites)*(2**Nsites),(cutoff**Nsites)*(2**Nsites)])).astype(np.complex64)

    # convert each term to matrix and sum up
    for t in expr.args:
        #timer=Timer(str(t)+' to matrix')
        #timer.start()
        termMat=convert_term_to_matrix(t,cutoff,Nsites, aops, adags, xs, xdags).astype(np.complex64)
        tmp=fullHam+termMat
        fullHam=tmp
        #timer.stop()
    # now need to drop all the buffers...
    #for i
        
    return fullHam


def convert_term_to_matrix(term, cutoff, Nsites, aops, adags, xs, xdags):
    #setupTimer=Timer('setup')
    #setupTimer.start()
    
    buffer=0
    prodMatrix = scipy.sparse.coo_matrix(np.eye(((cutoff+buffer)**Nsites)*(2**Nsites))).astype(np.complex64)
    
    #print("prod Start=",prodMatrix)
    
    if len(term.args)==0:
        return term*prodMatrix

    coef=1
    start=0
    
    if getattr(term.args[0],'__module__', None)=='sympy.core.numbers':
        coef=term.args[0]
        start=1

    #print(prodMatrix.shape)
    siteSubs = site_subs(cutoff+buffer, Nsites, aops, adags, xs, xdags)
    
    #setupTimer.stop()
    #print(term)
    
    # more efficient way is to group terms by site, multiple those small matrices
    # then kron product them.
    
    for t in term.args[start:]:
        #tSetupTimer=Timer(str(t) + 'setup')
        #tSetupTimer.start()
        isPow=False
        isBoson=False
        
        if type(t)==sp.core.power.Pow:
            isPow=True
        
        site=None
        siteMatrix=None
        
        if isPow:
            siteMatrix = siteSubs[t.args[0]]
            site = t.args[0].site
            siteMatrix = np.linalg.matrix_power(siteMatrix,t.args[1])
            if t.args[0].name[0]=='a':
                isBoson=True
        else:
            siteMatrix = siteSubs[t]
            site = t.site
            if t.name[0]=='a':
                isBoson=True
            
        #print("site={}, siteMatrix={}".format(site, siteMatrix))    
        #tSetupTimer.stop()
        
        siteMatrix=scipy.sparse.coo_matrix(siteMatrix).astype(np.complex64)
        
        #tFullTimer=Timer(str(t) + 'full')
        #tFullTimer.start()
        fullMatrix=scipy.sparse.coo_matrix([1]).astype(np.complex64)
        
        for i in range(0,Nsites):            
            fermId=scipy.sparse.coo_matrix(np.eye(2)).astype(np.complex64)
            bosId=scipy.sparse.coo_matrix(np.eye(cutoff+buffer)).astype(np.complex64)
            
            if site==str(i):
                if isBoson:
                    tmp=scipy.sparse.kron(fullMatrix,siteMatrix)
                    tmp2=scipy.sparse.kron(tmp,fermId)
                    fullMatrix=tmp2
                else:
                    tmp=scipy.sparse.kron(fullMatrix,bosId)
                    tmp2=scipy.sparse.kron(tmp,siteMatrix)
                    fullMatrix=tmp2
            else: 
                tmp=scipy.sparse.kron(fullMatrix,bosId)
                tmp2=scipy.sparse.kron(tmp,fermId)
                fullMatrix=tmp2
            
        fullMatrix=fullMatrix.astype(np.complex64)
        #tFullTimer.stop()
        
        
        
        ##print("fullMatrix={}".format(fullMatrix))
        #tMultTimer=Timer(str(t)+'multiply')
        #tMultTimer.start()
        
        prodMatrix=prodMatrix.dot(fullMatrix)
        
        #tMultTimer.stop()
        #print("prodMatrix={}".format(prodMatrix))
    
    
    # TODO if there's a buffer what's the right way 
    # to slice the matrix.
    return coef*prodMatrix


