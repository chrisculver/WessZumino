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
        allSubs[a[n]]=aMat(cutoff)
        allSubs[adag[n]]=aDagMat(cutoff)
        allSubs[x[n]]=xMat()
        allSubs[xd[n]]=xDagMat()
        
    return allSubs
    



def convert_to_matrix(expr, cutoff, buffer):
    
    new_expr = np.zeros([cutoff+buffer,cutoff+buffer])

    if type(expr)==sp.core.add.Add:
        for elem in expr.args:
            new_expr = new_expr + convert_term_to_matrix(elem,cutoff+buffer)
            
    elif type(expr)==sp.core.mul.Mul:
        tmp=convert_term_to_matrix(expr,cutoff+buffer)
        new_expr = new_expr+tmp
        
    elif type(expr)==sp.core.numbers.Float:
        new_expr = new_expr + convert_term_to_matrix(expr,cutoff+buffer)
    else:
        raise ValueError('Cannot convert type {} to matrix'.format(type(expr)))
        
        
    return np.array(new_expr.tolist())[:cutoff,:cutoff]





def convert_term_to_matrix(term, cutoff, Nsites):

    new_elem = np.eye(cutoff)
    has_aadag = False
    for elem in term.args:
        for i in sp.preorder_traversal(elem):
            if( type(i) is sp.tensor.indexed.Indexed ):
                has_aadag=True
    
    if has_aadag and type(term)==sp.core.mul.Mul:
        for elem in term.args:
            is_operator = False

            for i in sp.preorder_traversal(elem):
                if( type(i) is sp.tensor.indexed.Indexed ):
                    is_operator = True
                
            if is_operator:
                lst = elem.args

                print("ERROR: this needs to be done correctly...")
                
                if(len(lst)>1): # this should handle pow?
                    for i in range(lst[1]):
                        new_elem=np.matmul(new_elem,lst[0].subs(site_subs(cutoff, Nsites)))
                else:

                    new_elem=np.matmul(new_elem,elem.subs(site_subs(cutoff, Nsites))) 
            else:

                new_elem=elem*new_elem
        
        return new_elem
    
    elif has_aadag and type(term)!=sp.core.mul.Mul:
        raise ValueError('Havent implemented doing convert_to_matrix for type {}'.format(type(term)))
    
    else:
        if(term.args==()):
            new_elem=new_elem*term
        else:
            for elem in term.args:
                new_elem=new_elem*elem
        return new_elem*np.eye(cutoff)