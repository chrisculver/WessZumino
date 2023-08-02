#!/usr/bin/env python
# coding: utf-8

import math
import os
import time
import sys
sys.path.append('../..')

from src.wess_zumino_model import WessZuminoModel
from src.matrix_to_ps import matrix_to_pse
from src.binary_encodings import standard_encode
from src.qiskit_utilities import pauli_dict_to_op, op_to_trotter

from qiskit import Aer
from qiskit.utils import QuantumInstance
from qiskit.algorithms.eigensolvers import VQD
from qiskit.algorithms.optimizers import COBYLA
from qiskit.circuit.library import RealAmplitudes
from qiskit import BasicAer
from qiskit import qpy
from qiskit import QuantumCircuit, transpile
from qiskit.primitives import Estimator, Sampler
from qiskit.algorithms.state_fidelities import ComputeUncompute

import json
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments:
#   Choice of boundary conditions (dirichlet, periodic, ...)
#   Choice of prepotential (only linear and quadratic implemented so far)
#   Number of sites
#   Boson cutoff
#   Number of eigenvalues to find
#   Tolerance of VQD
#   Label for output pdf
#   TODO: maxiter; betas; ansatz reps; potential parameters...
if len(sys.argv) < 8:
    print("Usage: %s " % str(sys.argv[0]), end='')
    print("<BCs> <potential> <num_sites> <cutoff> ", end='')
    print("<num_eigvals> <tol> <count>")
    sys.exit(1)

BCs = str(sys.argv[1])
pot_tag = str(sys.argv[2])
N = int(sys.argv[3])
cutoff = int(sys.argv[4])
k = int(sys.argv[5])
maxiter=10000
target = float(sys.argv[6])
tag = round(math.log10(target))
run = int(sys.argv[7])
if k == 3:      # TODO: Generalize
    betas=[2,2,2]
elif k==5:
    betas=[2,2,2,2,2]
else:
    print("ERROR: Only set up for 3 and 5 eigenvalues so far, not %d" % k)
outfile = "%s/%s_N%d_L%d_k%d_tol%d-run%d.pdf" \
          % (BCs, pot_tag, N, cutoff, k, tag, run)

print("%s prepotential with %d sites, " % (pot_tag, N), end='')
print("%s BCs and cutoff %d" % (BCs, cutoff))
if k == 3:      # TODO: Generalize
    print("Hard-coded maxiter=10000 and betas=[2,2,2]")
elif k==5:
    print("Hard-coded maxiter=10000 and betas=[2,2,2,2,2]")
else:
    print("ERROR: Only set up for 3 and 5 eigenvalues so far, not %d" % k)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct hamiltonian as matrix, Pauli string and operator
def potential(self, n):
    if pot_tag == "linear":
        # m*q with m set below
        return self.mass*self.qs[n]
    elif pot_tag == "quad":
        # g*q^2 + c with g=1 and c=0
        return self.qs[n]*self.qs[n]
    else:
        print("ERROR: Unrecognized potential %s" % pot_tag)
        sys.exit(1)
    print("ERROR: Should never reach this part of 'potential' def")
    sys.exit(1)
    return np.nan * self.qs[n]

# Check BCs sanity
if BCs not in ['dirichlet', 'periodic']:
    print("ERROR: Unrecognized BCs %s" % BCs)
    sys.exit(1)

wz = WessZuminoModel(N, 1.0, potential, BCs)

runtime = -time.time()
wz.construct_ham_matrix(cutoff)
runtime += time.time()
print("%0.1f seconds to construct hamiltonian" % runtime)

runtime = -time.time()
ps=matrix_to_pse(wz.hamMat, standard_encode)
runtime += time.time()
print("%0.1f seconds to convert to Pauli string" % runtime)

runtime = -time.time()
op=pauli_dict_to_op(ps.to_dict())
runtime += time.time()
print("%0.1f seconds to convert to operator" % runtime)

nq=math.floor(math.log2(wz.hamMat.shape[0]))
if not math.log2(wz.hamMat.shape[0]).is_integer():
    nq+=1
print("nq = %d" % nq, flush=True)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Set up estimator and ansatz
estimator = Estimator()
sampler = Sampler()
fidelity = ComputeUncompute(sampler)

ansatz = RealAmplitudes(nq, entanglement='circular', reps=2)
#print(ansatz)

counts=[]
values=[]
steps=[]
def callback(eval_count, params, value, meta, step):
    counts.append(eval_count)
    values.append(value)
    steps.append(step)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Go
counts=[]
values=[]
steps=[]
runtime = -time.time()
vqd=VQD(estimator, fidelity, ansatz, optimizer=COBYLA(maxiter=maxiter, tol=target), k=k, betas=betas, callback=callback)
result=vqd.compute_eigenvalues(operator=op)
runtime += time.time()
print("%0.1f seconds for VQD" % runtime)

counts=np.asarray(counts)
steps=np.asarray(steps)
values=np.asarray(values)
for i in range(1,k+1):
    _counts=counts[np.where(steps==i)]
    _values=values[np.where(steps==i)]
    plt.plot(_counts,_values,label=r'$E_{{{}}}$'.format(i))
    if _counts[-1] == maxiter:
        print("NOT CONVERGED")
    print("E%d=%.8g after %d iterations" % (i, _values[-1], _counts[-1]))
plt.xlabel(r'Iteration')
plt.ylabel(r'$E$')
plt.yscale('log')
plt.legend(loc='upper right')
plt.savefig(outfile, bbox_inches='tight')
# ------------------------------------------------------------------
