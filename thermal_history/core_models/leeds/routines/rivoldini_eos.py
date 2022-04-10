"""
Routines provided by A. Rivoldini. Important function to use for model is liquidus (at bottom of file).

Description
------------------------------
The following code contains routines to compute the Equation of States (EoS)
for l-Fe-S and l-Fe-Si based on Terasaki 2019 (10.1029/2019JE005936).    
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy import integrate
from scipy import optimize


param = {'MFe':55.845,
         'MS':32.065,
         'MFeS':55.845+32.065,
         'RGas':8.31446}
     
class eosAndersonGrueneisen:
    def __init__(self,M0,p0,T0,V0,alpha0,KT0,KTP0,deltaT,kappa,
                 GibbsE=None,gamma0=None,q=None):
        self.pMax=200
        self.nbrPNodes=10001

        if (GibbsE is not None and (gamma0 is not None or q is not None)):
            print("Gibbs function and gamma not supported")
    				
        if (GibbsE is not None):
            self.GibbsFlag=True
            self.gamma0=0
            self.q=0
        else:
            self.GibbsFlag=False
            self.gamma0=gamma0
            self.q=q
    
        self.M0=M0
        self.p0=p0
        self.T0=T0
        self.V0=V0
        self.alpha0=alpha0
        self.KT0=KT0
        self.KTP0=KTP0
        self.deltaT=deltaT
        self.kappa=kappa
        self.GibbsE=GibbsE

        self.zetaA=np.zeros(self.nbrPNodes)
        self.px=np.zeros(self.nbrPNodes)
        self.zetaA[0]=1
        for i in range(1,self.nbrPNodes):
            self.px[i]=i/self.pMax
            self.zetaA[i]=self.compress(self.px[i])

        self.poly = CubicSpline(self.px,self.zetaA)

    def volume(self,x,T):
        # volume/V0
        p=x*self.pMax
        eta=(self.poly.__call__(p))**3
        alpha=self.alpha0*np.exp(-self.deltaT/self.kappa*(1-eta**self.kappa))
        return eta*np.exp(alpha*(T-self.T0))
    
    def Gibbs(self,p,T):
        if (p>self.p0):
            Gp = integrate.quad(lambda x: self.volume(x,T),
                                self.p0/self.pMax,p/self.pMax)[0]
        else :
            Gp=0
        return self.GibbsE(T)+1.e3*Gp*self.V0*self.pMax
        

    
    def compress(self,p):
        
        def VinetEq(x,p,KTP0,KT0):
            return -p+(3*np.exp((3*(-1+KTP0)*(1-x))/2)*KT0*(1-x))/x**2
            
        return optimize.brentq(VinetEq, 0.7, 1.2,args = (p,self.KTP0,self.KT0))
                                                                      

    def eos(self,p,T):
        deltaTemp=1 # temperature step for numerical differentiation, if too small results too noisy
        if (p>self.pMax):
            print("p should be smaller than ",self.pMax)
        T0=self.T0
        V0=self.V0
        alpha0=self.alpha0
        KT0=self.KT0
        KTP0=self.KTP0
        deltaT=self.deltaT
        kappa=self.kappa

        zeta=self.poly.__call__(p)
        eta=zeta**3
        alpha=alpha0*np.exp(-deltaT/kappa*(1-eta**kappa))
        V=V0*eta*np.exp(alpha*(T-T0))

        KT=(KT0*(4+(-5+3*KTP0)*zeta+3*(1-KTP0)*zeta**2))/np.exp((3*(-1+KTP0)*(-1+zeta))/2)
        KT=KT/(2*zeta**2)
        KT=KT/(1+(T-T0)*deltaT*alpha*eta**kappa)

        KTP=0.5*(KTP0-1)*zeta
        KTP=KTP+(8/3+(KTP0-5/3)*zeta)/(3*(4/3 +(KTP0-5/3)*zeta+(1-KTP0)*zeta**2))

        if (self.GibbsFlag):
            Gibbs=self.Gibbs(p,T)
            S=-(self.Gibbs(p,T+deltaTemp)-self.Gibbs(p,T-deltaTemp))/(2*deltaTemp)
            Cp=-T*(self.Gibbs(p,T+deltaTemp)-2*Gibbs+self.Gibbs(p,T-deltaTemp))/deltaTemp**2 # numerical second derivative of G with respect to T
            gamma=1/(Cp/(alpha*KT*V*1E+3)-alpha*T) # factor 1000 for conversion of GPa and cm^3/mol
            KS=KT*(1+gamma*alpha*T)
        else:
            Gibbs=0
            S=0
            gamma=self.gamma0*eta**self.q
            KS=KT*(1+gamma*alpha*T)
            Cp=1E+3*alpha*V*KS/gamma

        self.V=V
        self.rho=1.e3*self.M0/V
        self.alpha=alpha
        self.KT=KT
        self.KTP=KTP
        self.KS=KS
        self.gamma=gamma
        self.vp=np.sqrt(1E+9*KS/self.rho)
        self.vs=0
        self.Cp=Cp
        self.CV=1.e3*alpha*V*KT/gamma
        self.S=S
        self.GE=Gibbs	
        
        
    def __call__(self, p,T):
        self.eos(p,T)
        return {'M':self.M0,'V':self.V,'rho':self.rho,'alpha':self.alpha,'KT':self.KT,'KS': self.KS, 'gamma':self.gamma,'Cp': self.Cp, 'CV':self.CV,'vp':self.vp,'S':self.S,'GE':self.GE}        
        
               	     
class solution:
    def __init__(self,eos1,eos2,Vex):
        self.eM1=eos1
        self.eM2=eos2
        self.chi=np.zeros(2)
        self.Vex=Vex
        self.M0=0
        self.V=0
        self.Cp=0
        self.alpha=0
        self.KT=0
        self.gamma=0
        self.KS=0
        self.CV=0
        self.rho=0
        self.vp=0
        
    def __call__(self, x, p,T):
        self.eM1.eos(p,T)
        self.eM2.eos(p,T)
        chi = np.zeros(2)
        mFe=self.eM1.M0
        mS = self.eM2.M0-self.eM1.M0
        chi[1]=mFe*x/(mS+x*(mFe-mS)) # convert weight fraction to molar fraction S
        chi[1]=chi[1]/(1-chi[1]) # convert to molar fraction of FeS
        chi[0]=1-chi[1]
        Vexx=self.Vex(chi,p,T) #[Vex,dVex/dp,dVex/dT]
        M=[self.eM1.M0,self.eM2.M0]
        V=[self.eM1.V,self.eM2.V]
        Cp=[self.eM1.Cp,self.eM2.Cp]
        KT=[self.eM1.KT,self.eM2.KT]
        alpha=[self.eM1.alpha,self.eM2.alpha]
        self.M0=np.dot(M,chi)
        self.V=np.dot(V,chi)+Vexx[0]
        self.Cp=np.dot(Cp,chi)
        self.alpha=(np.dot([alpha[0]*V[0],alpha[1]*V[1]],chi)+Vexx[2])/self.V
        self.KT=-self.V/(-np.dot([V[0]/KT[0],V[1]/KT[1]],chi)+Vexx[1])
        self.gamma=1/(1E-3*self.Cp/(self.alpha*self.KT*self.V)-self.alpha*T)
        self.KS=self.KT*(1+self.alpha*self.gamma*T)
        self.CV=1E+3*self.alpha*self.V*self.KT/self.gamma
        self.rho=1E+3*self.M0/self.V
        self.vp=np.sqrt(1E+9*self.KS/self.rho)
        return {'M':self.M0,'V':self.V,'rho':self.rho,'alpha':self.alpha,'KT':self.KT,'KS': self.KS, 'gamma':self.gamma,'Cp': self.Cp, 'CV':self.CV,'vp':self.vp,'S':0,'GE':0}        
    
def TmFeA(p): # Fe liquidus from Anzellini et al. p in GPa
     return 495.4969600595926*(22.19 + p)**0.42016806722689076 

def TeFeS(p): # Fe-S eutectic to 13.5GPa
     return 1301.4062277227729 - 11.24327722772283*p

def xeFeS(p): # eutectic concetration 
    return 0.11 + 0.187*np.exp(-0.065*p)

def w2aS(x):
    MolarMassFe = param['MFe']
    MolarMassS = param['MS']
    return MolarMassFe*x/(MolarMassS+x*(MolarMassFe-MolarMassS))
    
class liquidusFeS:
    def __init__(self,deltaMu):
        self.deltaMu=deltaMu
        
    def __call__(self,x,p):
        return optimize.brentq(lambda T: self.deltaMu(w2aS(x),p,T), TeFeS(p)-20,TmFeA(p)+20) # bracket by eutectic and Fe liquidus
    

# Initialize the EoS 
fccFe=eosAndersonGrueneisen(M0=param["MFe"],p0=1.e-5,T0=298,V0=6.82,alpha0=7.e-5,KT0=163.4,KTP0=5.38,deltaT=5.5,kappa=1.4,GibbsE=(lambda T: 16300.921-395355.43/T-2476.28*np.sqrt(T)+ 381.47162*T+0.000177578*T**2-52.2754*T*np.log(T)))        
liquidFe=eosAndersonGrueneisen(M0=param['MFe'],p0=1E-5,T0=298,V0=6.88,alpha0=9E-5,KT0=148,KTP0=5.8,deltaT=5.1,kappa=0.56,GibbsE=(lambda T: 300-9007.3402+290.29866*T-46*T*np.log(T)))             
liquidFeS=eosAndersonGrueneisen(M0=param['MFeS'],p0=1E-5,T0=1650,V0=22.96,alpha0=11.9e-5, KT0=17.019,KTP0=5.922, deltaT=5.922,kappa=1.4,gamma0=1.3,q=0)

def VexFeFeS(chi,p,T):
    W11=-9.91275
    W12=0.731385
    W21=-1.32521
    W22=1.72716
    return (1-chi[1])*chi[1]*np.array([chi[1]*(W11+W12*np.log(1.5+p))+chi[0]*(W21+W22*np.log(1.5+p)),chi[1]*W12/(1.5+p)+chi[0]*W22/(1.5+p),0])
                                  
liquidCore=solution(liquidFe,liquidFeS,VexFeFeS) #Class to calculate properties from

def deltaMu(chi,p,T): # chemical potentals are equal at melting, assume no S in solid Fe
    y=chi/(1-chi) # mol fraction of Fe-S pseudo compound
    RGas=param['RGas']
    L1=53700.86652423554 - 3815.484280673008*p - 29.091372419282234*T
    L2=25339.70399255915 - 2951.106280178772*p
    fccFe.eos(p,T)
    liquidFe.eos(p,T)
    return -fccFe.GE+liquidFe.GE+RGas*T*np.log(1-y)+y**2 *(2*y-1)*L1 +2*(1-y)*y**2 *L2

#Melting temperature function to call from model. Takes mass fraction and pressure as arguments (see liquidusFeS).
liquidus=liquidusFeS(deltaMu)