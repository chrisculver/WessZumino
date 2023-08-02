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

import json
import numpy as np
import sympy as sp
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments:
#   Choice of boundary conditions (dirichlet, periodic, ...)
#   Choice of prepotential (only linear and quadratic implemented so far)
#   Number of sites
#   Boson cutoff
if len(sys.argv) < 5:
    print("Usage: %s " % str(sys.argv[0]), end='')
    print("<BCs> <potential> <num_sites> <cutoff>")
    sys.exit(1)

BCs = str(sys.argv[1])
pot_tag = str(sys.argv[2])
N = int(sys.argv[3])
cutoff = int(sys.argv[4])

print("%s prepotential with %d sites, " % (pot_tag, N), end='')
print("%s BCs and cutoff %d" % (BCs, cutoff))
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct hamiltonian as matrix, convert to Pauli strings
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

array = str(ps).split('+(') # TODO: Avoid '+' in complex coefficients
print("%d terms in Pauli string expression:" % len(array))
print("  %s" % array[0])
for i in range(1, len(array)):
  print("  (%s" % array[i])
# ------------------------------------------------------------------
