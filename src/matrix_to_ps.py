from sympy import symbols, expand, simplify
from collections import defaultdict
import scipy

def basis_to_pauli_string(i, j, q):
    if(i=='0' and j=='0'):
        return 0.5*symbols('I^'+q) + 0.5*symbols('Z^'+q)

    elif(i=='0' and j=='1'):
        return 0.5*symbols('X^'+q) + 0.5*1j*symbols('Y^'+q)

    elif(i=='1' and j=='0'):
        return 0.5*symbols('X^'+q) - 0.5*1j*symbols('Y^'+q)

    else:
        return 0.5*symbols('I^'+q) - 0.5*symbols('Z^'+q)


def matrix_to_pauli_strings(matrix, encoding):
    tst=scipy.sparse.coo_matrix(matrix)
    pauli_strings=0
    for i,j,v in zip(tst.row, tst.col, tst.data):
        pauli_strings += convert_element(v, i, j, matrix.shape[0], encoding)
        #pauli_strings=simplify(pauli_strings)
        # really slow to do for every element, maybe do every so often?
    return simplify(pauli_strings)


def convert_element(elem, i, j, N, encoding):
    nBinDigits = len(format(N-1,'b'))
    i_bin = encoding(i, nBinDigits)
    j_bin = encoding(j, nBinDigits)
    qubit_product = 1
       
    for q in range(0,len(i_bin)):
        qubit_product *= basis_to_pauli_string(i_bin[q], j_bin[q], str(len(i_bin)-1-q))
            
    return elem*qubit_product