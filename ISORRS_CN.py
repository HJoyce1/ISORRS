def Crank_Nicolson(dz,nz,dt,nt,D,U):
    import numpy as np
    from scipy.sparse import diags
    # dy - gridsize in space
    # ny - size of spatial grid
    # dt - gridsize for time
    # nt - number of time steps (always set to 2 in this circumstance)
    # D - coefficient (set to 1 in this circumstance)
    # U - input vector
    
    Uout = [] # list for storing V arrays at certain time steps
    U0 = U[0] # boundary condition on left side
    U1 = U[-1] # boundary condition on right side
    r = D*dt/dz**2  # diffusion number
    # create coefficient matrix:
    A = diags([-0.5*r, 1+r, -0.5*r], [-1, 0, 1], shape=(nz-2, nz-2)).toarray() 
    B1 = diags([0.5*r, 1-r, 0.5*r],[-1, 0, 1], shape=(nz-2, nz-2)).toarray()
    #breakpoint()
    #breakpoint()
    B = np.dot(U[1:-1],B1) #U[1:-1]
    B[0] = B[0]+0.5*r*(U0+U0)
    B[-1] = B[-1]+0.5*r*(U1+U1)
    Uout[1:-1] = np.linalg.solve(A,B)
    return Uout
