import sys
sys.path.append('..')

from setup import get_sympy_ham

from src.ham_to_sparse_matrix import *
from src.constants import *
from src.matrix_to_ps import matrix_to_pauli_strings
from src.binary_encodings import *
from src.timer import *
from src.qiskit_utilities import *

import sympy as sp
import scipy.sparse.linalg


N=2
hoHam,aops, adags, xs, xdags = get_sympy_ham(N)
cutoff=8 # 8 is more then a minute?
hamMat=convert_to_matrix(hoHam,cutoff,N,aops,adags,xs,xdags)

t=Timer("Matrix to PS")
t.start()
for i in range(10):
    ps=new_matrix_to_op(hamMat,standard_encode)
t.stop()
