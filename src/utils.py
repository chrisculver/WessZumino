def apply_bc(ham, bcType):

    # note this is exactly hardcoded for this finite difference method.
    boundaryConditions = {}
    if bcType == 'periodic':
        boundaryConditions = {q[-1]: q[N-1], q[N]: q[0],
                          x[N]: x[0], xd[N]: xd[0],
                         }
    elif bcType == 'dirichlet':
        boundaryConditions = {q[-1]: 0, q[N]: 0,
                          x[N]: 0, xd[N]: 0
                         }
        
    
    for i in range(0,N):
        ham+=totHam.subs(n,i).subs(boundaryConditions)
    # ham.subs(boundaryConditions).doit() # this doesn't work?
    
    return ham