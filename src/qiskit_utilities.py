from qiskit.opflow import I,X,Y,Z,Zero,PauliTrotterEvolution
from qiskit.quantum_info.operators import Pauli
from qiskit.circuit import Parameter

import scipy
from src.sympy_utilities import sympy_expr_to_list
from sympy import expand
from collections import defaultdict
from src.matrix_to_ps import basis_to_pauli_string,convert_element

char_to_op={'I': I, 'X': X, 'Y': Y, 'Z': Z}

def pauli_string_to_trotter_step(ps, time):
    return op_to_trotter(pauli_string_to_op(ps),time)

def pauli_string_to_op(ps):
    return pauli_dict_to_op(pauli_string_to_dict(ps))

def matrix_to_qiskit_op(matrix, encoding):
    tst=scipy.sparse.coo_matrix(matrix)
    pauli_strings=0
    for i,j,v in zip(tst.row, tst.col, tst.data):
        pauli_strings += pauli_string_to_op(convert_element(v, i, j, matrix.shape[0], encoding))
        pauli_strings
        #pauli_strings=simplify(pauli_strings)
        # really slow to do for every element, maybe do every so often?
    return pauli_strings.reduce()

def pauli_string_to_dict(ps):
    pauli_dict = defaultdict(complex)
    ps = expand(ps)
    #goes through all terms in pauli_strings (aka all a_i in sum_i a_i)
    for arg in ps.args:
        #print(arg)
        #each term is a product of a number * paulis, convert each elem to list
        arg_list=sympy_expr_to_list(arg)


        paulis = [ '' for i in range(0,len(arg_list)-1) ]
        start = 1
        coef = complex(arg_list[0])

        if(len(str(arg_list[1]))==1):
            coef*=1j
            start=2
            paulis.remove('')

        for ai in range(start, len(arg_list)):
            symbol = str(arg_list[ai])
            parts = symbol.split('^')
            paulis[ int(parts[1]) ] = parts[0]

        key = ''
        for i in paulis:
            key += i

        pauli_dict[key] += coef
        
    return pauli_dict


def pauli_dict_to_op(pd):
    expr = 0
    for pstring, coef in pd.items():
        term = char_to_op[pstring[0]]
        for char in pstring[1:]:
            term = term^char_to_op[char]
        expr += coef.real*term
    return expr

def op_to_trotter(op, time):
    evo_time = Parameter('t')
    evolution_op = (evo_time*op).exp_i()
    fixed_time = evolution_op.bind_parameters({evo_time: time})
    return PauliTrotterEvolution(trotter_mode='suzuki').convert(fixed_time)

def new_convert_element(elem, i, j, N, encoding):
    nBinDigits = len(format(N-1,'b'))
    i_bin = encoding(i, nBinDigits)
    j_bin = encoding(j, nBinDigits)
    qubit_product = new_basis_to_pauli_string(i_bin[len(i_bin)-1],j_bin[len(j_bin)-1])
    
    for q in range(len(i_bin)-2,-1,-1):
        qubit_product = qubit_product^new_basis_to_pauli_string(i_bin[q], j_bin[q])

    return (complex(elem)*qubit_product).reduce()  


def new_basis_to_pauli_string(i, j):
    if(i=='0' and j=='0'):
        return 0.5*I + 0.5*Z

    elif(i=='0' and j=='1'):
        return 0.5*X + 0.5*1j*Y

    elif(i=='1' and j=='0'):
        return 0.5*X - 0.5*1j*Y

    else:
        return 0.5*I - 0.5*Z


def new_matrix_to_op(matrix,encoding):
    tst=scipy.sparse.coo_matrix(matrix)
    pauli_strings=0
    for i,j,v in zip(tst.row, tst.col, tst.data):
        pauli_strings+=new_convert_element(v, i, j, matrix.shape[0], encoding)
        pauli_strings.reduce()

    return pauli_strings.reduce()