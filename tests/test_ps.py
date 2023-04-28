import sys 
sys.path.append('..')

from src.pauli_strings import I,X,Y,Z,PauliStringExpr,PauliString


a=PauliString(1.0,[I,X,X,I])
b=PauliString(-1.0,[I,X,X,I])

assert PauliStringExpr({b.gates: b.coef})+b==0