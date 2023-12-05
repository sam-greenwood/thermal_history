from numba import njit
import numpy as np
from ...leeds_thermal.routines.diffusion import thomas_algorithm, diffusion, LHS, RHS


def diffusion_uneven(y, x, dt, D, k, dk_dx, lower_bc, upper_bc, constant=0, coord='spherical'):

    '''Numerical diffusion solution on an unevenly spaced grid. D and k must both be either a single number or 
    both an array of same size as x/y. dk_dx is either a single number or a single number/array if D/k are arrays.

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
        either 'cart' to denote cartesian coordinate system or 'spherical' for spherical coordinates, by default 'spherical'

    Returns
    -------
    array 
        Diffused values after timestep dt.
    '''

    dx = np.zeros(x.size)
    
    dx[1:] = np.diff(x) #Node spacing x(i)-x(i-1)
    dx[0] = dx[1]

    if np.array(D).size == 1:
        #All must be single numbers
        D     = np.ones(x.size)*D
        k     = np.ones(x.size)*k
        dk_dx = np.zeros(x.size)*dk_dx

    if np.array(dk_dx).size == 1:
        #Expand to an array if a single number is supplied
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
    lower, main, upper = LHS_uneven(x, dx, dt, D, k, dk_dx, l_type, u_type, coord=coord)

    #Solve for n iterations
    y_new = y.copy()
    for _ in range(n):
        B = RHS_uneven(y_new, x, dt, D, k, dk_dx, dx, lower_bc, upper_bc, coord=coord, constant=constant)
        y_new = thomas_algorithm(B, lower, main, upper)

    return y_new


@njit
def LHS_uneven(x, dx, dt, D, k, dk_dx, l_type, u_type, coord='spherical'):
    '''Construct left hand side of linear equations

    Parameters
    ----------
    x : array
        position array
    dx : array
        array containing grid spacings.
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
def RHS_uneven(y, x, dt, D, k, dk_dx, dx, lower_bc, upper_bc, coord='spherical', constant=0):
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
    dx : array
        array containing grid spacings
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
        upper[i] +=  dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)

        #Spherical coordinates term
        if coord == 'spherical':

            # #If lower boundary is at r=0, use l’Hospital’s rule to deal with singularity
            if x[i] == 0:
                lower[0] = 0
                main[0]  = -3*dt*D[i]/(dx[i]**2) + 1
                upper[0] =  3*dt*D[i]/(dx[i]**2)
            else:
                lower[i] += -dt*D[i]/x[i] / (dx[i]+dx2)
                upper[i] +=  dt*D[i]/x[i] / (dx[i]+dx2)

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





def diffusion_discont(y, x, x_discont, dt,  D_lower, D_upper, k_lower, k_upper, lower_bc, upper_bc, constant=0, coord='spherical'):

    '''Numerical diffusion solution on an unevenly spaced grid for a discontinuous conductivity.

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
        either 'cart' to denote cartesian coordinate system or 'spherical' for spherical coordinates, by default 'spherical'

    Returns
    -------
    array 
        Diffused values after timestep dt.
    '''

    dx = np.zeros(x.size)
    
    dx[1:] = np.diff(x) #Node spacing x(i)-x(i-1)
    dx[0] = dx[1]


    m = (np.max([D_lower, D_upper])*dt/(2*np.min(dx)**2)) #uses minimum grid spacing present to estimate accurate time stepping.

    #Check solution will be stable. Run multiple times with smaller dt if necessary
    n=1
    if np.max(m) > 0.5:
        n = int(np.ceil(np.max(m)/0.5))
        dt = dt/n

    #Get boundary condition types from tuples of BC types and values
    l_type = lower_bc[0]
    u_type = upper_bc[0]

    #Discritise into linear equations (Ay=B) under Crank-Nicolson scheme.
    lower, main, upper = LHS_discont(x, x_discont, dx, dt, D_lower, D_upper, k_lower, k_upper, l_type, u_type, coord=coord)

    #Solve for n iterations
    y_new = y.copy()
    for _ in range(n):
        B = RHS_discont(y_new, x, x_discont, dt, D_lower, D_upper, k_lower, k_upper, dx, lower_bc, upper_bc, coord=coord, constant=constant)
        y_new = thomas_algorithm(B, lower, main, upper)

    return y_new



def LHS_discont(x, x_discont, dx, dt, D_lower, D_upper, k_lower, k_upper, l_type, u_type, coord='spherical'):
    '''Construct left hand side of linear equations

    Parameters
    ----------
    x : array
        position array
    dx : array
        array containing grid spacings.
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
        

        if x[i] == x_discont:

            # Difference in heat flow across interface
            lower[i] = k_lower/(x[i]-x[i-1])
            main[i]  = -(k_upper/(x[i+1]-x[i]) + k_lower/(x[i]-x[i-1]))
            upper[i] = k_upper/(x[i+1]-x[i])

        else:

            if x[i] < x_discont:
                D = D_lower
                k = k_lower
            else:
                D = D_upper
                k_upper


            #Standard diffusion term
            lower[i] = -dt*D/(dx[i]*(dx[i]+dx2))
            main[i]  = 1 + dt*D/(dx[i]*dx2)
            upper[i] = -dt*D/(dx2*(dx[i]+dx2))

            # No Variable conductivity yet
            # #Add terms for variable conductivity
            # lower[i] += dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)
            # upper[i] += -dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)

            #Spherical coordinates term
            if coord == 'spherical':

                # #If lower boundary is at r=0, use l’Hospital’s rule to deal with singularity
                if x[i] == 0:
                    lower[0] = 0
                    main[0]  = 3*dt*D/(dx[i]**2) + 1
                    upper[0] = -3*dt*D/(dx[i]**2)

                else:
                    lower[i] += dt*D/x[i] / (dx[i]+dx2)
                    upper[i] += -dt*D/x[i] / (dx[i]+dx2)

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



def RHS_discont(y, x, x_discont, dt, D_lower, D_upper, k_lower, k_upper, dx, lower_bc, upper_bc, coord='spherical', constant=0):
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
    dx : array
        array containing grid spacings
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


        if x[i] == x_discont:

            # No difference in heat flow across interface.
            lower[i] = 0
            main[i] = 0
            upper[i] = 0

        else:

            if x[i] < x_discont:
                D = D_lower
                k = k_lower
            else:
                D = D_upper
                k = k_upper
        
            #Standard diffusion term
            lower[i] = dt*D/(dx[i]*(dx[i]+dx2))
            main[i]  = 1 - dt*D/(dx[i]*dx2)
            upper[i] = dt*D/(dx2*(dx[i]+dx2))

            # No variable conductivity for now
            # #Add terms for variable conductivity
            # lower[i] += -dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)
            # upper[i] +=  dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)

            #Spherical coordinates term
            if coord == 'spherical':

                # #If lower boundary is at r=0, use l’Hospital’s rule to deal with singularity
                if x[i] == 0:
                    lower[0] = 0
                    main[0]  = -3*dt*D/(dx[i]**2) + 1
                    upper[0] =  3*dt*D/(dx[i]**2)
                else:
                    lower[i] += -dt*D/x[i] / (dx[i]+dx2)
                    upper[i] +=  dt*D/x[i] / (dx[i]+dx2)

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



def diffusion_discont_variable(y, x, x_discont, dt, D, k, dk_dx, lower_bc, upper_bc, constant=0, coord='spherical'):

    '''Numerical diffusion solution on an unevenly spaced grid, with a discontinuity in the conductivity. D and k must both be either a single number or 
    both an array of same size as x/y. dk_dx is either a single number or a single number/array if D/k are arrays.

    Parameters
    ----------
    y : array
        Values to diffuse
    x : array
        position grid
    x_discont : float
        x position of the discontinuity. Should be a node at this location.
    dt : float
        time step
    D : array
        diffusivity
    k : array
        thermal conductivity
    dk_dx : array
        spatial derivative of conductivity. Inf value at x_discont will be ignored.
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
        either 'cart' to denote cartesian coordinate system or 'spherical' for spherical coordinates, by default 'spherical'

    Returns
    -------
    array 
        Diffused values after timestep dt.
    '''

    dx = np.zeros(x.size)
    
    dx[1:] = np.diff(x) #Node spacing x(i)-x(i-1)
    dx[0] = dx[1]

    if np.array(D).size == 1:
        #All must be single numbers
        D     = np.ones(x.size)*D
        k     = np.ones(x.size)*k
        dk_dx = np.zeros(x.size)*dk_dx

    if np.array(dk_dx).size == 1:
        #Expand to an array if a single number is supplied
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
    lower, main, upper = LHS_discont_variable(x, dx, x_discont, dt, D, k, dk_dx, l_type, u_type, coord=coord)

    #Solve for n iterations
    y_new = y.copy()
    for _ in range(n):
        B = RHS_discont_variable(y_new, x, x_discont, dt, D, k, dk_dx, dx, lower_bc, upper_bc, coord=coord, constant=constant)
        y_new = thomas_algorithm(B, lower, main, upper)

    return y_new


@njit
def LHS_discont_variable(x, dx, x_discont, dt, D, k, dk_dx, l_type, u_type, coord='spherical'):
    '''Construct left hand side of linear equations

    Parameters
    ----------
    x : array
        position array
    dx : array
        array containing grid spacings.
    x_discont : float
        x position of the discontinuity
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
        
        if x[i] == x_discont:

            k_lower, k_upper = k[i-1], k[i+1]
            # Difference in heat flow across interface
            lower[i] = k_lower/(x[i]-x[i-1])
            main[i]  = -(k_upper/(x[i+1]-x[i]) + k_lower/(x[i]-x[i-1]))
            upper[i] = k_upper/(x[i+1]-x[i])
        
        else:

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
def RHS_discont_variable(y, x, x_discont, dt, D, k, dk_dx, dx, lower_bc, upper_bc, coord='spherical', constant=0):
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
    dx : array
        array containing grid spacings
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

        if x[i] == x_discont:

            # No difference in heat flow across interface.
            lower[i] = 0
            main[i] = 0
            upper[i] = 0
        
        else:

            #Standard diffusion term
            lower[i] = dt*D[i]/(dx[i]*(dx[i]+dx2))
            main[i]  = 1 - dt*D[i]/(dx[i]*dx2)
            upper[i] = dt*D[i]/(dx2*(dx[i]+dx2))

            #Add terms for variable conductivity
            lower[i] += -dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)
            upper[i] +=  dt*D[i]/(2*k[i]) * dk_dx[i] / (dx[i]+dx2)

            #Spherical coordinates term
            if coord == 'spherical':

                # #If lower boundary is at r=0, use l’Hospital’s rule to deal with singularity
                if x[i] == 0:
                    lower[0] = 0
                    main[0]  = -3*dt*D[i]/(dx[i]**2) + 1
                    upper[0] =  3*dt*D[i]/(dx[i]**2)
                else:
                    lower[i] += -dt*D[i]/x[i] / (dx[i]+dx2)
                    upper[i] +=  dt*D[i]/x[i] / (dx[i]+dx2)

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


