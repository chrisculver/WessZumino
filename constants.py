import sympy as sp

#pn, qn, qnP1, qnM1 = sp.symbols('p_n, q_n, q_{n+1}, q_{n-1}', commutative=False)
#xn, xdn, xnP1, xdnP1 = sp.symbols('\chi_n', '\chi^{\dagger}_n', '\chi_{n+1}', '\chi^{\dagger}_{n+1}')
n=sp.symbols('n')
aLat=sp.symbols('a')
V=sp.Function('V')
p=sp.IndexedBase('p',commutative=False)
q=sp.IndexedBase('q',commutative=False)
x=sp.IndexedBase('\chi',commutative=False)
xd=sp.IndexedBase('\chi^{\dagger}',commutative=False)
a = sp.IndexedBase('a',commutative=False)
adag = sp.IndexedBase('a^{\dagger}',commutative=False)