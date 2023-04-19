import sys
sys.path.append('..')
from src.ham_to_sparse_matrix import *
from src.constants import *
from src.matrix_to_ps import matrix_to_pauli_strings
from src.binary_encodings import *
from src.timer import *

import sympy as sp


def get_sympy_ham(N):
    bosonNI=sp.expand( pn**2/(2*aLat) + (aLat/2)*((qnP1-qnM1)/(2*aLat))**2 )
    bosonI=sp.expand( (aLat/2)*V(qn)**2 + aLat*V(qn)*(qnP1-qnM1)/(4*aLat) + aLat*(qnP1-qnM1)*V(qn)/(4*aLat) )
    fermionNI=sp.expand( -(xdnP1*xn+xdn*xnP1)/(2*aLat) )
    fermionI=sp.expand( sp.diff(V(qn),qn)*(xdn*xn-(1/2)) )

    aVal=1

    # depends on finite-difference method
    qs=[SiteSymbol('q',str(i)) for i in range(-1,N+1,1)]
    # maybe make qs a normal site list
    # and make an extra boundaryQs list for q[-1], q[N]
    # counting would be normal computer science way for rest of code.

    ps=[SiteSymbol('p',str(i)) for i in range(-1,N+1,1)] # don't really need extras
    aops=[SiteSymbol('a',str(i)) for i in range(-1,N+1,1)]
    adags=[SiteSymbol('a^{\dagger}',str(i)) for i in range(-1,N+1,1)]
    xs=[SiteSymbol('\chi',str(i)) for i in range(-1,N+1,1)]
    xdags=[SiteSymbol('\chi^{\dagger}',str(i)) for i in range(-1,N+1,1)]

    # note this is exactly hardcoded for this finite difference method.
    bcType = 'dirichlet'
    boundaryConditions = {}
    if bcType == 'periodic':
        boundaryConditions = {qs[0]: qs[N], qs[N+1]: qs[1],
                            xs[0]: -xs[N], xs[N+1]: -xs[1], 
                            xdags[0]: -xdags[N],  xdags[N+1]: -xdags[1]
                            }
        
    elif bcType == 'dirichlet':
        boundaryConditions = { qs[0]: 0, qs[N+1]: 0,
                            xs[0]: 0, xs[N+1]: 0,
                            xdags[0]: 0, xdags[N+1]:0
                            }

    m=1

    ham=0

    mass=1

    for i in range(1,N+1):
        ham+=(bosonNI+bosonI+fermionNI).subs({
            pn: ps[i],
            qn: qs[i], qnP1: qs[i+1], qnM1: qs[i-1],
            xn: xs[i], xnP1: xs[i+1],
            xdn: xdags[i], xdnP1: xdags[i+1]
        }).subs(boundaryConditions)

    def potential(n):
        # m*q with m=1
        return -mass*qs[n]

    potentialSubs={}
    for n in range(1,N+1):
        potentialSubs[V(qs[n])]=potential(n)


    ham=sp.simplify(ham.subs(potentialSubs).subs(aLat,aVal))

    for i in range(1,N+1):
        if i%2==0:
            ham+=sp.simplify(fermionI.subs({qn: qs[i], xn: xs[i], xdn: xdags[i]}).subs(potentialSubs))
        else:
            ham-=sp.simplify(fermionI.subs({qn: qs[i], xn: xs[i], xdn: xdags[i]}).subs(potentialSubs))
    ham=sp.simplify(ham.subs(aLat,aVal))
    sp.expand(ham)

    HOdofSubs = {}
    #offset because of BC
    for i in range(1,N+1):
        HOdofSubs[qs[i]] = 0.5*sp.sqrt(2/m)*(aops[i] + adags[i])
        HOdofSubs[ps[i]] = complex(0,1)*sp.sqrt(2*m)*(adags[i] - aops[i])/2 

    hoHam=sp.expand(ham.subs(HOdofSubs))
    hoHam=sp.nsimplify(hoHam,tolerance=1e-8)
    
    return hoHam, aops, adags, xs, xdags