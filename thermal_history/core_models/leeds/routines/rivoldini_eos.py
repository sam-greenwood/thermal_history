#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description
------------------------------
The following code contains routines to compute the Equation of States (EoS)
for l-Fe-S based on Terasaki 2019 (10.1029/2019JE005936).    
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy import integrate
from scipy import optimize


param = {'MFe':55.845,
         'MS':32.065,
         'MC':12.011,
         'MFeS':55.845+32.065,
         'RGas':8.31446}

class eosDinsdale:
    def __init__(self,M0,V0,KT0,KTP0,a0,a1,a,b,c,d,e,f,g):
        self.M0=M0
        self.V0=V0
        self.KT0=KT0*1.e9
        self.KTP0=KTP0
        self.a0=a0
        self.a1=a1
        self.a=a
        self.b=b
        self.c=c
        self.d=d
        self.e=e
        self.f=f
        self.g=g
    
    def eos(self,px,T):
        M0,V0,KT0,KTP0=self.M0,self.V0,self.KT0,self.KTP0
        a0,a1=self.a0,self.a1
        a,b,c,d,e,f,g=self.a,self.b,self.c,self.d,self.e,self.f,self.g
        p=1.e9*px
        self.GE=a + g/T**3 + f/T**2 +e/T+b*T+d*T**2 +(np.exp(a0*T+(a1*T**2)/2.)*KT0*(-1+(1+(KTP0*p)/KT0)**((-1 + KTP0)/KTP0))*V0)/(-1 + KTP0)+c*T*np.log(T)
        self.S=-2*d - (12*g)/T**5 - (6*f)/T**4 - (2*e)/T**3-c/T-(np.exp(a0*T+(a1*T**2)/2.)*KT0*(-1+(1+(KTP0*p)/KT0)**((-1+KTP0)/KTP0))*(a0**2+a1+2*a0*a1*T+a1**2*T**2)*V0)/(-1+KTP0)
        self.Cp=-(T*(2*d+(12*g)/T**5+(6*f)/T**4+(2*e)/T**3+c/T+(np.exp(a0*T+(a1*T**2)/2.)*KT0*(-1+(1+(KTP0*p)/KT0)**((-1+KTP0)/KTP0))*(a0**2+a1+2*a0*a1*T+a1**2*T**2)*V0)/(-1+KTP0)))
        self.V=((np.exp(a0*T+(a1*T**2)/2.)*V0)/(1+(KTP0*p)/KT0)**(1/KTP0))
        self.alpha=a0+a1*T
        self.KT=KT0+KTP0*p
        self.gamma=1/(self.Cp/(self.alpha*self.KT*self.V)-self.alpha*T)
        self.KS=self.KT*(1+self.gamma*self.alpha*T)
        self.CV=self.alpha*self.V*self.KT/self.gamma
        self.V=1.e6*self.V
        self.rho=1.e3*self.M0/self.V
        self.vp=np.sqrt(self.KS/self.rho)
        self.KT=1.e-9*self.KT
        self.KS=1.e-9*self.KS
        
        
    def __call__(self, p,T):
        self.eos(p,T)
        return {'M':self.M0,'V':self.V,'rho':self.rho,'alpha':self.alpha,'KT':self.KT,'KS': self.KS, 'gamma':self.gamma,'Cp': self.Cp, 'CV':self.CV,'vp':self.vp,'S':self.S,'GE':self.GE}

class eosAndersonGrueneisen:
    def __init__(self,M0,p0,T0,V0,alpha0,KT0,KTP0,deltaT,kappa,
                 GibbsE=None,gamma0=None,q=None):
        self.pMax=200
        self.nbrPNodes=10001

        if (GibbsE is not None and (gamma0 is not None or q is not None)):
            print("Gibbs function and gamma not supported")
    				
        if (GibbsE is not None):
            self.GibbsFlag=True
            self.VZFlag=False
            self.gamma0=0
            self.q=0  
        elif (GibbsE is None and gamma0 is None and q is None):
            self.GibbsFlag=False
            self.VZFlag=True
            self.gamma0=0
            self.q=0
        else:
            self.GibbsFlag=False
            self.VZFlag=False
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
        elif (self.VZFlag):
            Gibbs=0
            S=0
            gamma=(1/2*KTP-5/6+2/9*p/KT)/(1-4/3*p/KT)
            KS=KT*(1+gamma*alpha*T)
            Cp=1E+3*alpha*V*KS/gamma
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
    def __init__(self,eos1,eos2,ex,FeX=True):
        self.eM1=eos1
        self.eM2=eos2
        self.chi=np.zeros(2)
        self.ex=ex
        self.M0=0
        self.V=0
        self.Cp=0
        self.alpha=0
        self.alphaC=0
        self.KT=0
        self.gamma=0
        self.KS=0
        self.CV=0
        self.rho=0
        self.vp=0
        self.FeX=FeX
        
    def __call__(self, x, p,T): # x is wt fraction of light element X
        self.eM1.eos(p,T) # Fe
        self.eM2.eos(p,T) # FeX or X
        MFe=self.eM1.M0
        chi = np.zeros(2)
        if (self.FeX):
            MX=self.eM2.M0-self.eM1.M0 # molar mass of X
            chi[1]=x/(1-x)*MFe/MX # convert weight fraction of X to molar fraction of FeX
            dxdchi=MFe*MX/(MFe+MX*chi[1])**2 # derivative of weight fraction with respect to molar fraction of X
        else:
            MX=self.eM2.M0
            chi[1]=MFe*x/(MX+x*(MFe-MX)) # convert weight fraction of X to molar fraction of X
            dxdchi=(MX*MFe)/(MFe*(-1+chi[1]) - MX*chi[1])**2            
        chi[0]=1-chi[1]
        ex=self.ex(chi,p,T)
        M=[self.eM1.M0,self.eM2.M0]
        V=[self.eM1.V,self.eM2.V]
        Cp=[self.eM1.Cp,self.eM2.Cp]
        KT=[self.eM1.KT,self.eM2.KT]
        alpha=[self.eM1.alpha,self.eM2.alpha]
        self.M0=np.dot(M,chi)
        self.V=np.dot(V,chi)+ex['V']
        self.Cp=np.dot(Cp,chi)
        self.alpha=(np.dot([alpha[0]*V[0],alpha[1]*V[1]],chi)+ex['dVdT'])/self.V
        if (self.FeX):
            self.alphaC=-(MX/self.M0-(ex['dVdchi']+self.eM2.V-self.eM1.V)/self.V)/dxdchi
        else:
            self.alphaC=-((MX-MFe)/self.M0-(ex['dVdchi']+self.eM2.V-self.eM1.V)/self.V)/dxdchi
        self.KT=-self.V/(-np.dot([V[0]/KT[0],V[1]/KT[1]],chi)+ex["dVdp"])
        self.gamma=1/(1E-3*self.Cp/(self.alpha*self.KT*self.V)-self.alpha*T)
        self.KS=self.KT*(1+self.alpha*self.gamma*T)
        self.CV=1E+3*self.alpha*self.V*self.KT/self.gamma
        self.rho=1E+3*self.M0/self.V
        self.vp=np.sqrt(1E+9*self.KS/self.rho)
        return {'M':self.M0,'V':self.V,'rho':self.rho,'alpha':self.alpha,'alphaC':self.alphaC,'KT':self.KT,'KS': self.KS, 'gamma':self.gamma,'Cp': self.Cp, 'CV':self.CV,'vp':self.vp,'S':0,'GE':0}        
    
def TmFeA(p): # Fe liquidus from Anzellini et al. p in GPa
     return 495.4969600595926*(22.19 + p)**0.42016806722689076 

def TeFeS(p): # Fe-S eutectic to 13.5GPa
     return 1301.4062277227729 - 11.24327722772283*p
     
def TmFeS(p): # FeS melting T from Boehler 1992
          return 774.1353011127796*(9.938473089891751 + p)**0.2800438973467849

def xeFeS(p): # eutectic concetration 
    return 0.11 + 0.187*np.exp(-0.065*p)

def TeFeS(p): #eutectic temperatur for Fe-S for p<=10GPa
    return 1348.1385906216763 - 25.746183402852235*p - 0.37593441022701973*p**2 + 0.11225946350331356*p**3

def deltSFeS(T): #entropy difference between solid and liquif FeS, data is from NIST at 1bar in J/K/mol
    return -191.00037241868853 - 35.5021515055/T**2 - 0.085561502*T + 0.000024360150336031002*T**2 - 3.368803973666667e-14*T**3 + 33.276999999999994*np.log(T)

def w2aS(x):
    MolarMassFe = param['MFe']
    MolarMassS = param['MS']
    return MolarMassFe*x/(MolarMassS+x*(MolarMassFe-MolarMassS))

def a2wS(chi):
    MolarMassFe = param['MFe']
    MolarMassS = param['MS']
    return MolarMassS*chi/((1-chi)*MolarMassFe+chi*MolarMassS)
    
def w2aC(x):
    MolarMassFe = param['MFe']
    MolarMassC = param['MC']
    return MolarMassFe*x/(MolarMassC+x*(MolarMassFe-MolarMassC))

def a2wC(chi):
    MolarMassFe = param['MFe']
    MolarMassC = param['MC']
    return MolarMassC*chi/((1-chi)*MolarMassFe+chi*MolarMassC)
        
    
# melting T of FeC and all that along the graphite stability region
class pdFeC():
    def Tm(self,x,p):
        return -973.1949862972292 - 159.45023271197547*p - 2.649550966097275*p**2 + 64212.43996334035*x + 1271.4356827385757*p*x - 222434.25817139592*x**2
    def dTmdp(self,x,p):
        return -159.45023271197562 - 5.299101932194669*p + 1271.4356827385748*x
    def dTmdx(self,x,p):
        return 64212.43996334037 + 1271.4356827385748*p - 444868.5163427913*x
    def xg(self,p): # lowest C at p where graphite is stable
        return 0.0437052 + 0.00355333*p
    def Txg(self,p): # lowest T at p where graphite is stable
        return 1431.19 + 45.3298*p
    
    
class liquidusFeS:
    def __init__(self,deltaMu):
        self.deltaMu=deltaMu
        
    def __call__(self,x,p):
        if (x<xeFeS(p)):
            return optimize.brentq(lambda T: self.deltaMu(w2aS(x),p,T), TeFeS(p)-20,TmFeA(p)+100) # bracket by eutectic and Fe liquidus
        else:
            return TeFeS(p)+(TmFeS(p)-TeFeS(p))*(xeFeS(p)-x)/(xeFeS(p)-a2wS(0.5))
        
def kFeSPommier(x,T): #x in wt
    return (2.445*T)/(92.71381967698574 + 1189.0330473455936*x + 2404.678175564527*x**2)

def kFeSWagle(x,T): # x in wt
    return (2.445*T)/(75 + 127*x + 327*x**2)

def kFeCZhang(x,T): # x in wt
    return (2.445*T)/(75 + (2420.6*x)/(0.27401+x))       
    
    
# Initialize EoS 
fccFe=eosAndersonGrueneisen(M0=param["MFe"],p0=1.e-5,T0=298,V0=6.82,alpha0=7.e-5,KT0=163.4,KTP0=5.38,deltaT=5.5,kappa=1.4,GibbsE=(lambda T: 16300.921-395355.43/T-2476.28*np.sqrt(T)+ 381.47162*T+0.000177578*T**2-52.2754*T*np.log(T)))        
liquidFe=eosAndersonGrueneisen(M0=param['MFe'],p0=1E-5,T0=298,V0=6.88,alpha0=9E-5,KT0=148,KTP0=5.8,deltaT=5.1,kappa=0.56,GibbsE=(lambda T: 300-9007.3402+290.29866*T-46*T*np.log(T)))                                                                                           
liquidFeS=eosAndersonGrueneisen(M0=param['MFeS'],p0=1E-5,T0=1650,V0=22.96,alpha0=11.9e-5, KT0=17.019,KTP0=5.922, deltaT=5.922,kappa=1.4,gamma0=1.3,q=0)

FeSV=eosAndersonGrueneisen(M0=param['MFeS'],p0=1E-5,T0=1000,V0=18.042,alpha0=10.42e-5, KT0=54.3,KTP0=4, deltaT=2.0679,kappa=1.4)

graphite=eosDinsdale(M0=param['MC'],V0=5.259e-6,KT0=33.33,KTP0=12,a0=2.32e-5,a1=5.7e-9,a=-17368.441,b=170.73,c=-24.3,d=-0.4723e-3,e=2562600,f=-2.643e8,g=1.2e10)
diamond=eosDinsdale(M0=param['MC'],V0=3.412e-6,KT0=588.235,KTP0=5,a0=2.43e-6,a1=1e-8,a=-17368.441+1009,b=170.73+4.88,c=-24.3-0.01,d=-0.4723e-3,e=2562600+135400.,f=-2.643e8+33.05e5,g=1.2e10-9.e8)
liquidC=eosDinsdale(M0=param['MC'],V0=3.089e-6,KT0=42.034,KTP0=10.1,a0=2.32e-5,a1=5.7e-9,a=-17368.441+117369,b=170.73-24.63,c=-24.3,d=-0.4723e-3,e=2562600,f=-2.643e8,g=1.2e10)


def exFeFeS(chi,p,T):
    W11=-9.91275
    W12=0.731385
    W21=-1.32521
    W22=1.72716
    
    W11 = -12.418572079527602
    W21 = -2.3752788584923823
    W22 = 1.9348238292453082
    W12 = 1.0196634224150105
    V=np.log(1.5+p)
    dVdp=1/(1.5+p)
    return {'V':(1-chi[1])*chi[1]*(chi[1]*(W11+W12*V)+chi[0]*(W21+W22*V)),
            'dVdp':(1-chi[1])*chi[1]*(chi[1]*W12*dVdp+chi[0]*W22*dVdp),
            'dVdT':0,
            'dVdchi':-(chi[1]*(-2+3*chi[1])*(W11+W12*V))+(-1+chi[1])*(-1+3*chi[1])*(W21+W22*V)}
    
                              
liquidCore=solution(liquidFe,liquidFeS,exFeFeS)

liquidCoreIdeal=solution(liquidFe,liquidFeS,lambda chi,p,T:{'V':0,'dVdp':0,'dVdT':0,'dVdchi':0})

liquidCoreFeCideal=solution(liquidFe,liquidC,lambda chi,p,T:{'V':0,'dVdp':0,'dVdT':0,'dVdchi':0},FeX=False)    

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

