from numba import njit
import numpy as np


def diffusion2(y, x, dt, D, k, dk_dx, lower_bc, upper_bc, constant=0, coord='spherical'):
    '''Numerical diffusion solution

    Parameters
    ----------
    y : array
        Values to diffuse
    x : array
        position grid
    dt : float
        time step
    D : float or array
        diffusivity
    k : float or array
        thermal conductivity
    dk_dx : float or array
        spatial derivative of conductivity
    lower_bc : tuple
        tuple of 2 values to denote lower boundary condition. First value must be 0 or 1 to denote
        condition on the 0th or 1st spatial derivative of the solution. (dirichlet or nuemann type).
        Second value is the value of the boundary condition
    upper_bc : tuple
        tuple of 2 values to denote upper boundary condition. First value must be 0 or 1 to denote
        condition on the 0th or 1st spatial derivative of the solution. (dirichlet or nuemann type).
        Second value is the value of the boundary condition
    constant : float or array, optional
        Constant in the partial differential eqution, by default 0
    coord : str, optional
        either 'cart' to denote cartesian coordinate system of 'spherical' for spherical coordinates, by default 'spherical'

    Returns
    -------
    array 
        Diffused values after timestep dt.
    '''

    dx = np.zeros(x.size)
    
    dx[1:] = np.diff(x) #Node spacing x(i)-x(i-1)
    dx[0] = dx[1]

    if np.array(D).size == 1:
        D     = np.ones(x.size)*D
        k     = np.ones(x.size)*k
        dk_dx = np.zeros(x.size)*D

    if np.array(dk_dx).size == 1:
        dk_dx = np.ones(x.size)*dk_dx

    m = (D*dt/(2*np.min(dx)**2)) #uses minimum grid spacing present to estimate accurate time stepping.

    #Check solution will be stable. Run multiple times with smaller dt if necessary
    n=1
    if np.max(m) > 0.5:
        n = int(np.ceil(np.max(m)/0.5))
        dt = dt/n

    #Get boundary condition types from tuples of BC types and values
    l_type = lower_bc[0]
    u_type = upper_bc[0]

    #Discritise into linear equations (Ay=B) under Crank-Nicolson scheme.
    lower, main, upper = LHS2(x, dx, dt, D, k, dk_dx, l_type, u_type, coord=coord)

    y_new = y.copy()
    for _ in range(n):
        B = RHS2(y_new, x, dt, D, k, dk_dx, dx, lower_bc, upper_bc, coord=coord, constant=constant)
        y_new = thomas_algorithm(B, lower, main, upper)

    return y_new


@njit
def LHS2(x, dx, dt, D, k, dk_dx, l_type, u_type, coord='spherical'):
    '''Construct left hand side of linear equations

    Parameters
    ----------
    x : array
        position array
    m : array
        m constant defined in diffusion()
    dt : float
        time step
    D : array
        diffusivity
    k : array
        conductivity
    dk_dx : array
        spatial derivative of conductivity
    l_type : int
        Lower boundary condition type (0: dirichlet, 1: neumann)
    u_type : int    
        Upper boundary condition type (0: dirichlet, 1: neumann)
    coord : str, optional
        either 'cart' to denote cartesian coordinate system of 'spherical' for spherical coordinates, by default 'spherical'

    Returns
    -------
    (array , array, array)
        lower, central, upper diagonals for left hand side in linear equations (A matrix in Ax=b)
    '''

    #Set diaganols in matrix
    lower = np.zeros(x.size)
    main = np.zeros(x.size)
    upper = np.zeros(x.size)

    for i in range(x.size):

        if i == x.size-1:
            dx2 = dx[i]
        else:
            dx2 = dx[i+1]
        
        #Standard diffusion term
        lower[i] = -dt*D[i]/(dx[i]*(dx[i]+dx2))
        main[i]  = 1 + dt*D[i]/(dx[i]*dx2)
        upper[i] = -dt*D[i]/(dx2*(dx[i]+dx2))


        #Add terms for variable conductivity
        lower[i] += dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)
        upper[i] += -dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)

        #Spherical coordinates term
        if coord == 'spherical':

            # #If lower boundary is at r=0, use l’Hospital’s rule to deal with singularity
            if x[i] == 0:
                lower[0] = 0
                main[0]  = 3*dt*D[i]/(dx[i]**2) + 1
                upper[0] = -3*dt*D[i]/(dx[i]**2)

            else:
                lower[i] += dt*D[i]/x[i] / (dx[i]+dx2)
                upper[i] += -dt*D[i]/x[i] / (dx[i]+dx2)

    #Handle boundary conditions
    if l_type == 0:
        upper[0] = 0
        main[0] = 1
    elif l_type == 1:
        upper[0] += lower[0]

    if u_type == 0:
        lower[-1] = 0
        main[-1] = 1
    elif u_type == 1:
        lower[-1] += upper[-1]

    return lower[1:], main, upper[:-1]


@njit
def RHS2(y, x, dt, D, k, dk_dx, dx, lower_bc, upper_bc, coord='spherical', constant=0):
    '''Construct right hand side of linear equations

    Parameters
    ----------
    y : array
        Values to diffuse
    x : array
        position array
    dt : array
        timestep
    D : array
        diffusivity
    k : array
        conductivity
    dk_dx : array
        spatial derivative of conductivity
    m : array
        m constant defined in diffusion()
    lower_bc : float
        Lower boundary condition value
    upper_bc : float    
        Upper boundary condition value
    coord : str, optional
        either 'cart' to denote cartesian coordinate system of 'spherical' for spherical coordinates, by default 'spherical'
    constant : float or array, optional
        Constant in the partial differential eqution, by default 0

    Returns
    -------
    array 
        Right hand side in linear equations (b vector in Ax=b)
    '''

    #unpack boundary conditions
    l_type, l_bc = lower_bc
    u_type, u_bc = upper_bc


    #Set diaganols in matrix
    lower = np.zeros(x.size)
    main = np.zeros(x.size)
    upper = np.zeros(x.size)


    for i in range(x.size):

        if i == x.size-1:
            dx2 = dx[i]
        else:
            dx2 = dx[i+1]
        
        #Standard diffusion term
        lower[i] = dt*D[i]/(dx[i]*(dx[i]+dx2))
        main[i]  = 1 - dt*D[i]/(dx[i]*dx2)
        upper[i] = dt*D[i]/(dx2*(dx[i]+dx2))

        #Add terms for variable conductivity
        lower[i] += -dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)
        upper[i] += dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)

        #Spherical coordinates term
        if coord == 'spherical':

            # #If lower boundary is at r=0, use l’Hospital’s rule to deal with singularity
            if x[i] == 0:
                lower[0] = 0
                main[0]  = -3*dt*D[i]/(dx[i]**2) + 1
                upper[0] = 3*dt*D[i]/(dx[i]**2)
            else:
                lower[i] += -dt*D[i]/x[i] / (dx[i]+dx2)
                upper[i] += dt*D[i]/x[i] / (dx[i]+dx2)

    C = np.diag(lower[1:],-1) + np.diag(main) + np.diag(upper[:-1],1)
    B = np.dot(C,y) + constant #NOTE use of constant has not been tested yet

    #Handle boundary conditions
    if l_type == 0:
        B[0] = l_bc

    elif l_type == 1:
        upper[0] += lower[0]
        B[0] = main[0]*y[0] + upper[0]*y[1] - 4*dx[0]*l_bc*lower[0]

    if u_type == 0:
        B[-1] = u_bc

    elif u_type == 1:
        lower[-1] += upper[-1]
        B[-1] = lower[-1]*y[-2] + main[-1]*y[-1] + 4*dx[-1]*u_bc*upper[-1]

    return B








############################################################################################


def diffusion(y, x, dt, D, k, dk_dx, lower_bc, upper_bc, constant=0, coord='spherical'):
    '''Numerical diffusion solution

    Parameters
    ----------
    y : array
        Values to diffuse
    x : array
        position grid
    dt : float
        time step
    D : float or array
        diffusivity
    k : float or array
        thermal conductivity
    dk_dx : float or array
        spatial derivative of conductivity
    lower_bc : tuple
        tuple of 2 values to denote lower boundary condition. First value must be 0 or 1 to denote
        condition on the 0th or 1st spatial derivative of the solution. (dirichlet or nuemann type).
        Second value is the value of the boundary condition
    upper_bc : tuple
        tuple of 2 values to denote upper boundary condition. First value must be 0 or 1 to denote
        condition on the 0th or 1st spatial derivative of the solution. (dirichlet or nuemann type).
        Second value is the value of the boundary condition
    constant : float or array, optional
        Constant in the partial differential eqution, by default 0
    coord : str, optional
        either 'cart' to denote cartesian coordinate system of 'spherical' for spherical coordinates, by default 'spherical'

    Returns
    -------
    array 
        Diffused values after timestep dt.
    '''

    dx = x[1]-x[0]

    if np.array(D).size == 1:
        D     = np.ones(x.size)*D
        k     = np.ones(x.size)*k
        dk_dx = np.zeros(x.size)*D

    if np.array(dk_dx).size == 1:
        dk_dx = np.ones(x.size)*dk_dx

    m = (D*dt/(2*dx**2))

    #Check solution will be stable. Run multiple times with smaller dt if necessary
    n=1
    if np.max(m) > 0.5:
        n = int(np.ceil(np.max(m)/0.5))
        dt = dt/n
        m = m/n

    #Get boundary condition types from tuples of BC types and values
    l_type = lower_bc[0]
    u_type = upper_bc[0]

    #Discritise into linear equations (Ay=B) under Crank-Nicolson scheme.
    lower, main, upper = LHS(x, m, D, k, dk_dx, l_type, u_type, coord=coord)

    y_new = y.copy()
    for _ in range(n):

        B = RHS(y_new, x, D, k, dk_dx, m, lower_bc, upper_bc, coord=coord, constant=constant)

        y_new = thomas_algorithm(B, lower, main, upper)
        # thomas_algorithm(B, lower.copy(), main.copy(), upper.copy(), lower.size, y_new)

    return y_new


@njit
def LHS(x, m, D, k, dk_dx, l_type, u_type, coord='spherical'):
    '''Construct left hand side of linear equations

    Parameters
    ----------
    x : array
        position array
    m : array
        m constant defined in diffusion()
    D : array
        diffusivity
    k : array
        conductivity
    dk_dx : array
        spatial derivative of conductivity
    l_type : int
        Lower boundary condition type (0: dirichlet, 1: neumann)
    u_type : int    
        Upper boundary condition type (0: dirichlet, 1: neumann)
    coord : str, optional
        either 'cart' to denote cartesian coordinate system of 'spherical' for spherical coordinates, by default 'spherical'

    Returns
    -------
    (array , array, array)
        lower, central, upper diagonals for left hand side in linear equations (A matrix in Ax=b)
    '''

    dx = x[1]-x[0]

    #Set diaganols in matrix
    lower = m*-1
    main  = 1+2*m
    upper = m*-1

    #add terms for variable thermal conductivity
    lower +=  m*dk_dx*dx/(2*k)
    upper += -m*dk_dx*dx/(2*k)


    #add terms for spherical coordinates if necessary
    if coord == 'spherical':

        #If lower boundary is at r=0, use l’Hospital’s rule to deal with singularity
        if x[0] == 0:
            main[0]  = (1+6*m[0])
            upper[0] = -6*m[0]
            lower[0] = 0

            lower[1:] +=  m[1:]*dx/x[1:]
            upper[1:] += -m[1:]*dx/x[1:]
        else:
            lower +=  m*dx/x
            upper += -m*dx/x


    #Handle boundary conditions
    if l_type == 0:
        upper[0] = 0
        main[0] = 1
    elif l_type == 1:
        upper[0] += lower[0]

    if u_type == 0:
        lower[-1] = 0
        main[-1] = 1
    elif u_type == 1:
        lower[-1] += upper[-1]

    return lower[1:], main, upper[:-1]


@njit
def RHS(y, x, D, k, dk_dx, m, lower_bc, upper_bc, coord='spherical', constant=0):
    '''Construct right hand side of linear equations

    Parameters
    ----------
    y : array
        Values to diffuse
    x : array
        position array
    D : array
        diffusivity
    k : array
        conductivity
    dk_dx : array
        spatial derivative of conductivity
    m : array
        m constant defined in diffusion()
    lower_bc : float
        Lower boundary condition value
    upper_bc : float    
        Upper boundary condition value
    coord : str, optional
        either 'cart' to denote cartesian coordinate system of 'spherical' for spherical coordinates, by default 'spherical'
    constant : float or array, optional
        Constant in the partial differential eqution, by default 0

    Returns
    -------
    array 
        Right hand side in linear equations (b vector in Ax=b)
    '''

    #unpack boundary conditions
    l_type, l_bc = lower_bc
    u_type, u_bc = upper_bc


    dx = x[1]-x[0]

    #Set diaganols in matrix
    lower = m*1
    main  = 1-2*m
    upper = m*1

    #add terms for variable diffusivity
    lower += -m*dk_dx*dx/(2*k)
    upper +=  m*dk_dx*dx/(2*k)

    #add terms for spherical coordinates if necessary
    if coord == 'spherical':

        #If lower boundary is at r=0, use l’Hospital’s rule to deal with singularity
        if x[0] == 0:
            main[0]  = (1-6*m[0])
            upper[0] = 6*m[0]
            lower[0] = 0

            lower[1:] += -m[1:]*dx/x[1:]
            upper[1:] +=  m[1:]*dx/x[1:]
        else:
            lower += -m*dx/x
            upper +=  m*dx/x

    C = np.diag(lower[1:],-1) + np.diag(main) + np.diag(upper[:-1],1)
    B = np.dot(C,y) + constant

    #Handle boundary conditions
    if l_type == 0:
        B[0] = l_bc

    elif l_type == 1:
        upper[0] += lower[0]
        B[0] = main[0]*y[0] + upper[0]*y[1] - 4*dx*l_bc*lower[0]


    if u_type == 0:
        B[-1] = u_bc

    elif u_type == 1:
        lower[-1] += upper[-1]
        B[-1] = lower[-1]*y[-2] + main[-1]*y[-1] + 4*dx*u_bc*upper[-1]

    return B



@njit
def thomas_algorithm(RHS, lower, main, upper):
    '''Jit compiled implementation of the thomas algorithm for solving
    linear equations.

    Parameters
    ----------
    RHS : array
        B vector in Ax=B
    lower : array
        lower diagonal values in A matrix
    main : array
        central diagonal values in A matrix
    upper : array
        upper diagonal values in A matrix

    Returns
    -------
    array
        Diffused values
    '''

    a, b, c, d = np.copy(lower), np.copy(main), np.copy(upper), np.copy(RHS)

    n = RHS.size

    for i in range(1, n):
        w = a[i-1]/b[i-1]
        b[i] = b[i] - w*c[i-1]
        d[i] = d[i] - w*d[i-1]

    y= np.zeros(n)
    y[-1] = d[-1]/b[-1]

    for i in range(n-2, -1, -1):
        y[i] = (d[i]-c[i]*y[i+1])/b[i]

    return y
