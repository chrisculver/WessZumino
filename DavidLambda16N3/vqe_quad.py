import sys, os
sys.path.append('..')

import sympy as sp
from src.ham_to_sparse_matrix import *
from src.constants import *
from src.matrix_to_ps import matrix_to_pauli_strings
from src.binary_encodings import *
import scipy.sparse.linalg

from src.qiskit_utilities import *
from qiskit import Aer
from qiskit.utils import QuantumInstance

from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import COBYLA
from qiskit.circuit.library import RealAmplitudes

from qiskit import qpy
from src.timer import *

import pickle

from qiskit import QuantumCircuit, transpile


# potential parameters 
c=-0.2
cFileString="cm0p2"   # reads as c=-0.2
c2=1 # DO NOT CHANGE

# Lambda - bosonic cutoff
cutoff=4

# how many VQE runs to do - prints out best result
vqeShots=10

# sites and lattice spacing
N=2
aVal=1
mass=1
m=1







backend = Aer.get_backend('statevector_simulator')
qinstance = QuantumInstance(backend, seed_simulator=2, seed_transpiler=2)

bosonNI=sp.expand( pn**2/(2*aLat) + (aLat/2)*((qnP1-qnM1)/(2*aLat))**2 )
bosonI=sp.expand( (aLat/2)*V(qn)**2 + aLat*V(qn)*(qnP1-qnM1)/(4*aLat) + aLat*(qnP1-qnM1)*V(qn)/(4*aLat) )
fermionNI=sp.expand( -(xdnP1*xn+xdn*xnP1)/(2*aLat) )
fermionI=sp.expand( sp.diff(V(qn),qn)*(xdn*xn-(1/2)) )


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


ham=0

for i in range(1,N+1):
    ham+=(bosonNI+bosonI+fermionNI).subs({
        pn: ps[i],
        qn: qs[i], qnP1: qs[i+1], qnM1: qs[i-1],
        xn: xs[i], xnP1: xs[i+1],
        xdn: xdags[i], xdnP1: xdags[i+1]
    }).subs(boundaryConditions)

def potential(n):
    return c + c2*qs[n]*qs[n]

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
hoHam



# DAVID - MAIN PART - Gets H as matrix, gets qiskit op, runs vqe


print("Lambda={}  |  N={}  |  c={}".format(cutoff, N, c))
hamMat=convert_to_matrix(hoHam,cutoff,N,aops,adags,xs,xdags)
ens=scipy.sparse.linalg.eigs(hamMat,k=6,sigma=0.0)[0]

opFileName="Data/op_quad_{}_lambda{}.pickle".format(cFileString,cutoff)

op=None
if os.path.isfile(opFileName):
    print("reading op file")
    with open(opFileName, "rb") as f:
        op = pickle.load(f)

else:
    print("computing op")
    ps=matrix_to_pauli_strings(hamMat,standard_encode)
    op = pauli_string_to_op(ps)

    with open(opFileName, "wb") as f:
        pickle.dump(op, f)


nq=math.floor(math.log2(hamMat.shape[0]))
if not math.log2(hamMat.shape[0]).is_integer():
    nq+=1

ansatz = RealAmplitudes(nq, reps=cutoff)


vqe = VQE(ansatz=ansatz, optimizer=COBYLA(), quantum_instance=qinstance)
results=[]
for i in range(vqeShots):
    results.append(vqe.compute_minimum_eigenvalue(op).eigenvalue.real)

print("{} & {:.2e} & {:.2e} \\\\".format(cutoff, np.min(ens).real, np.array(results).min()))

qc = QuantumCircuit(nq,nq)
qc.append(op_to_trotter(op,0.1), [i for i in range(nq)])
tmp = transpile(qc, basis_gates = ['cx', 'u1', 'u2', 'u3', 'H', 'X', 'Y', 'Z'])
print("gates for one trotter step = ", tmp.count_ops().get('cx'))