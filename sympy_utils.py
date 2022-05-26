from collections import defaultdict

import sympy as sp
import numpy as np

import copy
import math

pauliX = np.array([[0,1],[1,0]])
pauliY = np.array([[0,-1j],[1j,0]])
pauliZ = np.array([[1,0],[0,-1]])

a = sp.IndexedBase('a',commutative=False)
adag = sp.IndexedBase('a^{\dagger}',commutative=False)
x=sp.IndexedBase('\chi',commutative=False)
xd=sp.IndexedBase('\chi^{\dagger}',commutative=False)

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



def site_subs(cutoff, Nsites):
    allSubs = {}
    
    for n in range(0,Nsites): 
        allSubs[a[n]]=sp.Matrix(aMat(cutoff))
        allSubs[adag[n]]=sp.Matrix(aDagMat(cutoff))
        allSubs[x[n]]=sp.Matrix(xMat())
        allSubs[xd[n]]=sp.Matrix(xDagMat())
        
    return allSubs
    



def convert_to_matrix(expr, cutoff, Nsites):
    
    fullHam = np.zeros([(cutoff**Nsites)*(2**Nsites),(cutoff**Nsites)*(2**Nsites)]).astype(np.complex64)
    print(fullHam.shape)

    for t in expr.args:
        fullHam=fullHam+convert_term_to_matrix(t,cutoff,Nsites)
        
        
    return fullHam





def convert_term_to_matrix(term, cutoff, Nsites):
    coef=1
    if type(term.args[0])==sp.core.numbers.Float:
        coef=term.args[0]
    else:
        raise TypeError("Expected coef found type {}".format(type(term.args[0])))

    mats={}
    for t in term.args[1:]:
        boson=False
        isPow=False
        if type(t)==sp.core.power.Pow:
            isPow=True
            if hasattr(t.args[0],'name'):
                if t.args[0].name[0]=='a':
                    boson=True
        if hasattr(t,'name'):
            if t.name[0]=='a':
                boson=True
    
        siteSubs = site_subs(cutoff,Nsites)
    
        if boson:
            if isPow:
                mats['b'+str(t.args[0].indices[0])]=np.linalg.matrix_power(np.array(
                    t.subs(siteSubs).doit()).astype(np.complex64),t.args[1])
            else:
                mats['b'+str(t.indices[0])]=np.array(
                    t.subs(siteSubs).doit()).astype(np.complex64)
        else:
            if isPow:
                mats['f'+str(t.args[0].indices[0])]=np.linalg.matrix_power(np.array(
                t.subs(siteSubs).doit()).astype(np.complex64),t.args[1])
            else:
                mats['f'+str(t.indices[0])]=np.array(
                t.subs(siteSubs).doit()).astype(np.complex64)

    fullMat=1
    for i in range(Nsites):
        if 'b'+str(i) in mats:
            fullMat=np.kron(fullMat,mats['b'+str(i)])
        else:
            fullMat=np.kron(fullMat,np.eye(cutoff))
    
        if 'f'+str(i) in mats:
            fullMat=np.kron(fullMat,mats['f'+str(i)])
        else:
            fullMat=np.kron(fullMat,np.eye(2))

    return coef*fullMat    