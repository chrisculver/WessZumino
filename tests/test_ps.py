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

print(1.0*X*X + 2.0*X*Y)

print(1.0*X*X + 2.0*X*Y - 2.0*X*Y - 3.0*Z*Z + 4.0*X*X - (-1.0*X*X + 1.0*Z*Z))
print("\n---Testing mul---")
print( (1.0*X + 1.0*Z)*(1.0*Y + 1.0*Z) )
print( (1.0*Y + 1.0*Z)*(1.0*X + 1.0*Z))

