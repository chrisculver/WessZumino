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



def site_subs(cutoff, Nsites, aops, adags):
    allSubs = {}
    
    for n in range(-1,Nsites+1): 
        allSubs[aops[n]]=sp.Matrix(aMat(cutoff))
        allSubs[adags[n]]=sp.Matrix(aDagMat(cutoff))
        #allSubs[x[n]]=sp.Matrix(xMat())
        #allSubs[xd[n]]=sp.Matrix(xDagMat())
        
    return allSubs
    



def convert_to_matrix(expr, cutoff, Nsites, aops, adags):
    # start with a matrix of zeros
    fullHam = np.zeros([(cutoff**Nsites)*(2**Nsites),(cutoff**Nsites)*(2**Nsites)]).astype(np.complex64)

    # convert each term to matrix and sum up
    for t in expr.args:
        fullHam=fullHam+convert_term_to_matrix(t,cutoff,Nsites, aops, adags)
        
    # now need to drop all the buffers...
    #for i
        
    return fullHam





def convert_term_to_matrix(term, cutoff, Nsites, aops, adags):
    coef=1
    start=0
    
    if getattr(term.args[0],'__module__', None)=='sympy.core.numbers':
        coef=term.args[0]
        start=1
        
    mats={}
    for t in term.args[start:]:
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
        
        siteSubs = site_subs(cutoff,Nsites, aops, adags)
    
        if boson:
            if isPow:              
                siteSubs = site_subs(cutoff+t.args[1],Nsites, aops, adags)
                #print(siteSubs)
                mats['b'+str(t.args[0].site)]=np.linalg.matrix_power(np.array(
                        t.args[0].subs(siteSubs).doit()).astype(np.complex64),t.args[1])[:cutoff,:cutoff]
            
            else:
                siteSubs = site_subs(cutoff,Nsites, aops, adags)
                mats['b'+str(t.site)]=np.array(
                    t.subs(siteSubs).doit()).astype(np.complex64)
        else:
            if isPow:
                siteSubs = site_subs(cutoff,Nsites, aops, adags)
                print("Warning: raising grassman to a power")
                mats['f'+str(t.args[0].site)]=np.linalg.matrix_power(np.array(
                t.subs(siteSubs).doit()).astype(np.complex64),t.args[1])
            
            else:
                siteSubs = site_subs(cutoff,Nsites, aops, adags)
                mats['f'+str(t.site)]=np.array(
                t.subs(siteSubs).doit()).astype(np.complex64)

    #fullMat=np.eye((cutoff**Nsites)*(2**Nsites)).astype(np.complex64)
    fullMat=1
    for i in range(Nsites):
    #    if i==0:
    #        offset=0
    #    else:
    #        offset=((cutoff**i)*(2**i))
        
        if 'b'+str(i) in mats:
    #        fullMat[offset:offset+cutoff,offset:offset+cutoff]=mats['b'+str(i)]
            fullMat=np.kron(fullMat,mats['b'+str(i)])
        else:
        #    fullMat[offset:offset+cutoff,offset:offset+cutoff]=mats['b'+str(i)]
            fullMat=np.kron(fullMat,np.eye(cutoff))
    
        if 'f'+str(i) in mats:
        #    offset+=2
        #    fullMat[offset:offset+2,offset:offset+2]=mats['f'+str(i)]
            fullMat=np.kron(fullMat,mats['f'+str(i)])
        else:
            #fullMat=np.kron(fullMat,np.eye(2))
            fullMat=np.kron(fullMat,np.eye(2))

            
    return coef*fullMat    





def convert_boson_to_matrix(expr, cutoff, Nsites, aops, adags):
    # start with a matrix of zeros
    fullHam = np.zeros([(cutoff**Nsites),(cutoff**Nsites)]).astype(np.complex64)

    # convert each term to matrix and sum up
    for t in expr.args:
        fullHam=fullHam+convert_boson_term_to_matrix(t,cutoff,Nsites, aops, adags)
        
    # now need to drop all the buffers...
    #for i
        
    return fullHam





# This is WRONG!
# need to multiply from left to right no just kron everything.
# this doesn't handle ad0 a1 a0 a1d ad0^3 a1 a1d
def convert_boson_term_to_matrix(term, cutoff, Nsites, aops, adags):
    coef=1
    start=0
    
    if getattr(term.args[0],'__module__', None)=='sympy.core.numbers':
        coef=term.args[0]
        start=1
        
    mats={}
    for t in term.args[start:]:
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
        
        siteSubs = site_subs(cutoff,Nsites, aops, adags)
    
        if boson:
            if isPow:              
                siteSubs = site_subs(cutoff+t.args[1],Nsites, aops, adags)
                #print(siteSubs)
                mats['b'+str(t.args[0].site)]=np.linalg.matrix_power(np.array(
                        t.args[0].subs(siteSubs).doit()).astype(np.complex64),t.args[1])[:cutoff,:cutoff]
            
            else:
                siteSubs = site_subs(cutoff,Nsites, aops, adags)
                mats['b'+str(t.site)]=np.array(
                    t.subs(siteSubs).doit()).astype(np.complex64)
        else:
            print("non-boson term??!?!?!")
            
    #fullMat=np.eye((cutoff**Nsites)*(2**Nsites)).astype(np.complex64)
    fullMat=1
    for i in range(Nsites):
    #    if i==0:
    #        offset=0
    #    else:
    #        offset=((cutoff**i)*(2**i))
        if 'b'+str(i) in mats:
            print(mats['b'+str(i)])
    #        fullMat[offset:offset+cutoff,offset:offset+cutoff]=mats['b'+str(i)]
            fullMat=np.kron(fullMat,mats['b'+str(i)])
        else:
        #    fullMat[offset:offset+cutoff,offset:offset+cutoff]=mats['b'+str(i)]
            fullMat=np.kron(fullMat,np.eye(cutoff))
            
    return coef*fullMat    