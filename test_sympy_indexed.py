import sympy as sp

t=sp.IndexedBase('t',commutative=False)
u=sp.IndexedBase('u',commutative=False)

ts=[]
for i in range(3):
  ts.append(sp.Symbol('t_{}'.format(i),commutative=False))

a=sp.Symbol('a',commutative=False)
b=sp.Symbol('b',commutative=False)

print('t commutativity is {}'.format(t._assumptions['commutative']))
print('t[0] commutativity is {}'.format(t[0]._assumptions['commutative']))
print('u commutativity is {}'.format(u._assumptions['commutative']))
print('u[0] commutativity is {}'.format(u[0]._assumptions['commutative']))

#print(t[0]*u[0]+u[0]*t[0])
print(ts[0]*ts[1]+ts[1]*ts[0])

#t_0*u_0
print(a*b+b*a)
