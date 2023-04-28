import sys 
sys.path.append('..')

from src.pauli_strings import I,X,Y,Z,PauliStringExpr,PauliString


# Ideally I would write code like 
# expr = 1.5*X*Y*Y*Z + 1.5*I*Y*Y*Z 
# expr += 0.5*I*I*I*X + 2.5*I*X*X*I

# need _mul_ for PG*PG and float*PG
# need add for PExpr Pexpr and PS PS and PExpr PS
# always store as string?
print(1.0*X)
print(I*X)
print(2.0*X*I*Z*Y*X*X)

