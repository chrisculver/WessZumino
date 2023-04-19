import sys
sys.path.append('..')

from setup import get_sympy_ham

from src.ham_to_sparse_matrix import *
from src.constants import *
from src.matrix_to_ps import matrix_to_pauli_strings
from src.binary_encodings import *
from src.timer import *

import sympy as sp
import scipy.sparse.linalg


N=2
hoHam,aops, adags, xs, xdags = get_sympy_ham(N)
cutoff=4 # 8 is more then a minute?
hamMat=convert_to_matrix(hoHam,cutoff,N,aops,adags,xs,xdags)

t=Timer("Matrix to PS")
t.start()
ps=matrix_to_pauli_strings(hamMat,standard_encode)
t.stop()