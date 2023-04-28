from src.ham_to_sparse_matrix import convert_to_matrix
from src.constants import *

import sympy as sp

class WessZuminoModel:
    # sets up the lattice Wess-Zumino model and constructs 
    # the explicit hamiltonian for $N$ sites with the given potential.
    def __init__(self, N, mass, potential, bcType):
        self.N=N
        self.mass=mass
        self.potential=potential
        self.bcType=bcType

        bosonNI=sp.expand( pn**2/(2*aLat) + (aLat/2)*((qnP1-qnM1)/(2*aLat))**2 )
        bosonI=sp.expand( (aLat/2)*V(qn)**2 + aLat*V(qn)*(qnP1-qnM1)/(4*aLat) + aLat*(qnP1-qnM1)*V(qn)/(4*aLat) )
        fermionNI=sp.expand( -(xdnP1*xn+xdn*xnP1)/(2*aLat) )
        fermionI=sp.expand( sp.diff(V(qn),qn)*(xdn*xn-(1/2)) )

        self.aVal=1  # value of lattice spacing

        # depends on finite-difference method
        # store qs for the potential
        self.qs=[SiteSymbol('q',str(i)) for i in range(-1,N+1,1)]
        # maybe make qs a normal site list
        # and make an extra boundaryQs list for q[-1], q[N]
        # counting would be normal computer science way for rest of code.
        ps=[SiteSymbol('p',str(i)) for i in range(-1,N+1,1)] # don't really need extras
        # ps and qs can be discarded since we work in harmonic oscillator basis.

        # we need to store all these symbols to construct the Hamiltonian as a matrix
        self.xs=[SiteSymbol('\chi',str(i)) for i in range(-1,N+1,1)]
        self.xdags=[SiteSymbol('\chi^{\dagger}',str(i)) for i in range(-1,N+1,1)]

        # note this is exactly hardcoded for this finite difference method.
        # different finite difference methods may require more rules
        boundaryConditions = {}
        if self.bcType == 'periodic':
            boundaryConditions = {self.qs[0]: self.qs[N], self.qs[N+1]: self.qs[1],
                                self.xs[0]: -self.xs[N], self.xs[N+1]: -self.xs[1], 
                                self.xdags[0]: -self.xdags[N],  self.xdags[N+1]: -self.xdags[1]
                                }
            
        elif self.bcType == 'dirichlet':
            boundaryConditions = { self.qs[0]: 0, self.qs[N+1]: 0,
                                self.xs[0]: 0, self.xs[N+1]: 0,
                                self.xdags[0]: 0, self.xdags[N+1]:0
                                }

        self.ham = 0

        for i in range(1,N+1):
            self.ham+=(bosonNI+bosonI+fermionNI).subs({
                pn: ps[i],
                qn: self.qs[i], qnP1: self.qs[i+1], qnM1: self.qs[i-1],
                xn: self.xs[i], xnP1: self.xs[i+1],
                xdn: self.xdags[i], xdnP1: self.xdags[i+1]
            }).subs(boundaryConditions)

        potentialSubs={}
        for n in range(1,N+1):
            potentialSubs[V(self.qs[n])]=self.potential(self,n)

        self.ham = sp.simplify(self.ham.subs(potentialSubs).subs(aLat, self.aVal))

        for i in range(1,N+1):
            if i%2==0:
                self.ham+=sp.simplify(fermionI.subs({qn: self.qs[i], xn: self.xs[i], xdn: self.xdags[i]}).subs(potentialSubs))
            else:
                self.ham-=sp.simplify(fermionI.subs({qn: self.qs[i], xn: self.xs[i], xdn: self.xdags[i]}).subs(potentialSubs))
        self.ham=sp.simplify(self.ham.subs(aLat,self.aVal))
        sp.expand(self.ham)

        # again we need to store aops and dags for matrix construction.
        self.aops=[SiteSymbol('a',str(i)) for i in range(-1,N+1,1)]
        self.adags=[SiteSymbol('a^{\dagger}',str(i)) for i in range(-1,N+1,1)]
        HOdofSubs = {}
        #offset because of BC
        for i in range(1,N+1):
            HOdofSubs[self.qs[i]] = 0.5*sp.sqrt(2/self.mass)*(self.aops[i] + self.adags[i])
            HOdofSubs[ps[i]] = complex(0,1)*sp.sqrt(2*self.mass)*(self.adags[i] - self.aops[i])/2 

        self.hoHam=sp.expand(self.ham.subs(HOdofSubs))
        self.hoHam=sp.nsimplify(self.hoHam,tolerance=1e-8)

    def construct_ham_matrix(self,cutoff):
        self.cutoff=cutoff
        self.hamMat=convert_to_matrix(self.hoHam, cutoff, self.N, self.aops, self.adags, self.xs, self.xdags)