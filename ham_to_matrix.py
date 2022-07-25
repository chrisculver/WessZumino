from collections import defaultdict

import sympy as sp
import numpy as np

import copy
import math

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
    
    for n in range(-1,Nsites+1): 
        allSubs[aops[n]]=sp.Matrix(aMat(cutoff))
        allSubs[adags[n]]=sp.Matrix(aDagMat(cutoff))
        allSubs[xs[n]]=sp.Matrix(xMat())
        allSubs[xdags[n]]=sp.Matrix(xDagMat())
        
    return allSubs
    

def convert_to_matrix(expr, cutoff, Nsites, aops, adags, xs, xdags):
    # start with a matrix of zeros
    fullHam = np.zeros([(cutoff**Nsites)*(2**Nsites),(cutoff**Nsites)*(2**Nsites)]).astype(np.complex64)

    # convert each term to matrix and sum up
    for t in expr.args:
        fullHam=fullHam+convert_term_to_matrix(t,cutoff,Nsites, aops, adags, xs, xdags)
        
    # now need to drop all the buffers...
    #for i
        
    return fullHam


def convert_term_to_matrix(term, cutoff, Nsites, aops, adags, xs, xdags):
    coef=1
    start=0
    
    if getattr(term.args[0],'__module__', None)=='sympy.core.numbers':
        coef=term.args[0]
        start=1

    buffer=0
    prodMatrix = np.eye(((cutoff+buffer)**Nsites)*(2**Nsites)).astype(np.complex64)
    print(prodMatrix.shape)
    siteSubs = site_subs(cutoff+buffer, Nsites, aops, adags, xs, xdags)
    
    #print(term)
    
    # more efficient way is to group terms by site, multiple those small matrices
    # then kron product them.
    for t in term.args[start:]:
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
        
        fullMatrix=1
        for i in range(0,Nsites):            
            if site==str(i):
                if isBoson:
                    fullMatrix=np.kron(fullMatrix,siteMatrix)
                    fullMatrix=np.kron(fullMatrix,np.eye(2))
                else:
                    fullMatrix=np.kron(np.eye(cutoff+buffer))
                    fullMatrix=np.kron(fullMatrix,siteMatrix)
            else:
                fullMatrix=np.kron(fullMatrix,np.eye(cutoff+buffer))
                fullMatrix=np.kron(fullMatrix,np.eye(2))
            
        
        #print("fullMatrix={}".format(fullMatrix))

        prodMatrix=np.matmul(prodMatrix,fullMatrix)
        
        #print("prodMatrix={}".format(prodMatrix))
    
    
    # TODO if there's a buffer what's the right way 
    # to slice the matrix.
    return coef*prodMatrix






def convert_boson_to_matrix(expr, cutoff, Nsites, aops, adags):
    # start with a matrix of zeros
    fullHam = np.zeros([(cutoff**Nsites),(cutoff**Nsites)]).astype(np.complex64)

    # convert each term to matrix and sum up
    for t in expr.args:
        fullHam=fullHam+convert_boson_term_to_matrix(t,cutoff,Nsites, aops, adags)
        
    # now need to drop all the buffers...
    #for i
        
    return fullHam





# This is inefficient.
def convert_boson_term_to_matrix(term, cutoff, Nsites, aops, adags):
    coef=1
    start=0
    
    if getattr(term.args[0],'__module__', None)=='sympy.core.numbers':
        coef=term.args[0]
        start=1

    buffer=0
    prodMatrix = np.eye((cutoff+buffer)**Nsites).astype(np.complex64)
    siteSubs = site_subs(cutoff+buffer, Nsites, aops, adags)
    
    #print(term)
    for t in term.args[start:]:
        
        isPow=False
        if type(t)==sp.core.power.Pow:
            isPow=True
        
        site=None
        siteMatrix=None
        
        if isPow:
            siteMatrix = siteSubs[t.args[0]]
            site = t.args[0].site
            siteMatrix = np.linalg.matrix_power(siteMatrix,t.args[1])
        else:
            siteMatrix = siteSubs[t]
            site = t.site
            
        #print("site={}, siteMatrix={}".format(site, siteMatrix))    
        
        fullMatrix=1
        for i in range(0,Nsites):
            if site==str(i):
                fullMatrix=np.kron(fullMatrix,siteMatrix)
            else:
                fullMatrix=np.kron(fullMatrix,np.eye(cutoff+buffer))
        
        #print("fullMatrix={}".format(fullMatrix))
        
        prodMatrix=np.matmul(prodMatrix,fullMatrix)
        
        #print("prodMatrix={}".format(prodMatrix))
    
    
    # TODO if there's a buffer what's the right way 
    # to slice the matrix.
    return coef*prodMatrix
            
            