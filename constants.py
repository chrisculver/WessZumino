import sympy as sp

#pn, qn, qnP1, qnM1 = sp.symbols('p_n, q_n, q_{n+1}, q_{n-1}', commutative=False)
#xn, xdn, xnP1, xdnP1 = sp.symbols('\chi_n', '\chi^{\dagger}_n', '\chi_{n+1}', '\chi^{\dagger}_{n+1}')



        
n=sp.symbols('n')
aLat=sp.symbols('a')
V=sp.Function('V')

# my custom class with description attribute
class SiteSymbol(sp.Symbol):
    def __new__(cls, name, site):
        obj = sp.Symbol.__new__(cls, name+'_{'+str(site)+'}',commutative=False)
        obj.site=site
        return obj


# make new objects with description
pn = SiteSymbol('p','n')
qn = SiteSymbol('q','n')
qnP1 = SiteSymbol('q','n+1')
qnM1 = SiteSymbol('q','n-1')

xn=SiteSymbol('\chi','n')
xdn=SiteSymbol('\chi^{\dagger}','n')
xdnP1=SiteSymbol('\chi^{\dagger}','n+1')
