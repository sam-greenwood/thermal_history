'''
Stand alone script to pre-compute values for melting temperature of Fe-S alloys using Rivoldini's EOS, then fitting with polynomials
in order to have a fast calculation of Tm for P=4-10GPa and Sulphur content (S) of 0-15 wt%.

Polynomial coefficients are defined in the following way:

Tm0 = x[0]      + x[1]*S   + x[2]*S**2 + .... x[N]*S**N
Tm1 = x[N+1]    + x[N+2]*S + ...
TmN = x[N(N+1)] + ....

Tm = Tm0 + Tm1*P + ... TmN*P**N

The coefficients of zero order in S (x[0], x[N+1]...) are defined by a fit to the pure iron parameterisation of Anzellini
Tm_fe = 495.4969600595926*(22.19 + P)**0.42016806722689076

The max polynomial degree N is set by default to 3 but can be changed if required.

The coefficients are calculated and printed to screen and also saved in a file called 'optimal_solution.txt'.

A plot is produced of the error between the polynomial representation and the EOS values of TM is plotted at the end by randomly sampling the P/S space.
'''

N = 3 #Max degree of fit

import numpy as np
from scipy.optimize import curve_fit
import sys

import rivoldini_eos as eos

P_array = np.linspace(4,10, 30)*1e9 #Range of Pressures to consider
S_array = np.linspace(0,0.15, P_array.size) #Range of Sulphur concentrations to consider
Tm = np.zeros((P_array.size, S_array.size))

# Initialize the EoS 
fccFe=eos.eosAndersonGrueneisen(M0=eos.param["MFe"],p0=1.e-5,T0=298,V0=6.82,alpha0=7.e-5,KT0=163.4,KTP0=5.38,deltaT=5.5,kappa=1.4,GibbsE=(lambda T: 16300.921-395355.43/T-2476.28*np.sqrt(T)+ 381.47162*T+0.000177578*T**2-52.2754*T*np.log(T)))        
liquidFe=eos.eosAndersonGrueneisen(M0=eos.param['MFe'],p0=1E-5,T0=298,V0=6.88,alpha0=9E-5,KT0=148,KTP0=5.8,deltaT=5.1,kappa=0.56,GibbsE=(lambda T: 300-9007.3402+290.29866*T-46*T*np.log(T)))                                                                                           

def deltaMu(chi,p,T): # chemical potentals are equal at melting, assume no S in solid Fe
    y=chi/(1-chi) # mol fraction of Fe-S pseudo compound
    RGas=eos.param['RGas']
    L1=53700.86652423554 - 3815.484280673008*p - 29.091372419282234*T
    L2=25339.70399255915 - 2951.106280178772*p
    fccFe.eos(p,T)
    liquidFe.eos(p,T)
    return -fccFe.GE+liquidFe.GE+RGas*T*np.log(1-y)+y**2 *(2*y-1)*L1 +2*(1-y)*y**2 *L2

liquidus=eos.liquidusFeS(deltaMu)



#Fit pure Fe melting curve of Anzellini. These will be coefficients when Sulphur concentration=0
x_fe = np.polyfit(P_array, 495.4969600595926*(22.19 + P_array/1e9)**0.42016806722689076, N)[::-1]

def forward_fit(u,*args):
    #Forward problem. Calcluate Tm using polynomial representation of Tm for
    #give S content and pressure. Takes coefficients for pure iron above
    #for the non S dependent terms.

    #Tm0 = x_fe[0] + x[1]*S + x[2]*S**2 + .... 
    #Tm1 = ....
    #TmN = ....

    #Tm = Tm0 + Tm1*P + ... TmN*P**N

    x = args

    S, P = u #Sulphur comp and Pressure in first input

    Tm = np.zeros(S.size)
    c = 0
    for i in range(N+1):

        temp = np.dot(x_fe[i], S**0) #Pure iron coefficient

        for j in range(1,N+1):
            temp += np.dot(x[c+j-1], S**j)

        Tm += temp*P**i
        c += N

    return Tm


def forward(u,*args):
    #Forward problem including pure iron coefficients in given parameters.

    x = args

    S, P = u #Sulphur comp and Pressure in first input

    Tm = np.zeros(S.size)
    c = 0
    for i in range(N+1):

        temp = np.zeros(S.size)

        for j in range(0,N+1):
            temp += np.dot(x[c+j], S**j)

        Tm += temp*P**i
        c += N+1

    return Tm


#Check to see if -f flag given to reuse saved coefficients in file.
flag = True
if len(sys.argv) > 1:
    if sys.argv[1] == '-f':
        flag = False
    else:
        print('only accepted arugment is -f, assuming no flags were given')

if flag:
    #Calculate Tm from eos at each P/S combination
    for i,P in enumerate(P_array):

        for j,S in enumerate(S_array):

            eutectic = eos.xeFeS(P/1e9)
            assert S <= eutectic, f'S concentration ({S}) above eutectic ({eutectic})' #Double check we stay on iron side of eutectic comp.

            if S == 0:
                Tm[i,j] = eos.TmFeA(P/1e9) #Pure iron melting temp
            else:
                Tm[i,j] = liquidus(S,P/1e9) #Alloy melting temp
        

        print(f'computing Tm at pressure {P/1e9} GPa')


    #Mesh data for passing into curve_fit
    S_mesh, P_mesh = np.meshgrid(S_array, P_array)
    xdata = np.vstack((S_mesh.ravel(), P_mesh.ravel()))

    x_initial = np.ones((N+1)**2-(N+1)) #Initial guess of ones for every coefficient

    #Fit 2D data to produce coefficients.
    xopt, xcov = curve_fit(forward_fit, xdata, Tm.ravel(), x_initial)

    xopt = np.insert(xopt, 0, x_fe[0])
    for i in range(1,N+1):
        xopt = np.insert(xopt, i*(N+1), x_fe[i])

    #Save solution to file to reuse
    with open('optimal_solution.txt', 'w') as f:
        for x in xopt:
            f.write(f'{x}\n')

else:
    xopt = np.genfromtxt('optimal_solution.txt') #Read in solution from file

print('Optimal Solution')
print(xopt)
print('------------------')

#Error across random points

S_rand = np.random.rand(150)*S_array[-1]
P_rand = np.random.rand(S_rand.size)*(P_array[-1]-P_array[0]) + P_array[0]
Tm_rand = np.zeros(S_rand.size)
Tm_predict = np.zeros(S_rand.size)


#Calculate Tm from eos and polynomial fit at random points in S/P space.
for i in range(S_rand.size):
    S, P = S_rand[i], P_rand[i]

    if S == 0:
        Tm_rand[i] = eos.TmFeA(P/1e9) #Pure iron melting temp
    else:
        Tm_rand[i] = liquidus(S,P/1e9) #Alloy melting temp

    Tm_predict[i] = forward((S,P), *xopt)

try:
    import matplotlib.pyplot as plt

    #Plot random points coloured by error
    y = Tm_predict-Tm_rand
    plt.scatter(S_rand*100, P_rand, s=200, c=y, marker='o', cmap='seismic', edgecolor='k', vmin=-np.max(np.abs(y)), vmax=np.max(np.abs(y)))
    plt.colorbar()
    plt.show()

    #Plot error between eos and polynomial as a function of S
    plt.plot(S_array*100, forward((S_array,P_array[0]*np.ones(S_array.size)), *xopt) - np.array([liquidus(S,P_array[0]/1e9) for S in S_array]))
    plt.xlabel('Sulphur composition [wt%]')
    plt.ylabel('Error in poly representation to eos.')
    plt.show()
except:
    pass