#!/usr/bin/env python
# coding: utf-8

import math
import os
import time
import sys
sys.path.append('../..')

from src.wess_zumino_model import WessZuminoModel

import json
import numpy as np
import sympy as sp
import scipy.linalg
import scipy.sparse.linalg
from scipy.sparse import csr_matrix
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments:
#   Choice of boundary conditions (dirichlet, periodic, ...)
#   Choice of prepotential (only linear and quadratic implemented so far)
#   Number of sites
#   Boson cutoff
#   Number of eigenvalues to find
if len(sys.argv) < 6:
    print("Usage: %s " % str(sys.argv[0]), end='')
    print("<BCs> <potential> <num_sites> <cutoff> <num_eigvals>")
    sys.exit(1)

BCs = str(sys.argv[1])
pot_tag = str(sys.argv[2])
N = int(sys.argv[3])
cutoff = int(sys.argv[4])
k = int(sys.argv[5])

print("%s prepotential with %d sites, " % (pot_tag, N), end='')
print("%s BCs and cutoff %d" % (BCs, cutoff))
print("Sparse computation of lowest %d energies" % k)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct hamiltonian as matrix
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
temp = scipy.sparse.linalg.eigs(wz.hamMat, k=k, sigma=0.0)[0]
runtime += time.time()
print("%0.1f seconds for sparse solve" % runtime)
eigs = np.sort(temp)
for i in range(k):
    print("  Sparse e%d = %.8g" % (i, eigs[i].real))
    if np.abs(eigs[i].imag) > 1e-4:
      print("  Warning: Im(e%d) = %.4g" % (i, eigs[i].imag))
#print("Turned off sparse solve to avoid singular factor runtime error")

runtime = -time.time()
dense = csr_matrix(wz.hamMat)
temp = scipy.linalg.eigvals(dense.todense())
runtime += time.time()
print("%0.1f seconds for full solve" % runtime)
eigs = np.sort(temp)
for i in range(k):
    print("  Full e%d = %.8g" % (i, eigs[i].real))
    if np.abs(eigs[i].imag) > 1e-4:
      print("  Warning: Im(e%d) = %.4g" % (i, eigs[i].imag))
#print("Turned off full solve to avoid overwhelming system")
# ------------------------------------------------------------------
