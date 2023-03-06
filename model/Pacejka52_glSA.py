#Pacejka Modell 52 Paramter
import numpy as np

#Parameter importieren
import parameter.Pacejka_Params_Indy as pa

#i increase for accuracy 
i=1000
SL=np.linspace(-1,1,i)
alpha = np.linspace(-15,15,i)
Fz = 8000
Fz0 = 14000
IA=0.1
#Fz0= PacLat_Fznom #nominal load (z.b 80% von Fzmax)

#Fx longitudenal
def Pacejka52_long(SL,Fz,Fz0,PHX1, PHX2, PKX1, PKX2, PKX3, PCX1, PDX1, PDX2, PDX3, PEX1, PEX2, PEX4, PVX1, PVX2, LMUX, LEX): 
    #IA=inclination angle,  Fz0= highest load in data, 
    IA = 0
    #converts to rad
    gamma = IA*np.pi/180
    
    dfz = ((Fz-Fz0)/Fz0)  
    kappa = SL #- Reifengeschw/ Fahrzeuggeschw
    SHx =  PHX1 +  PHX2 * dfz
    kappax = kappa + SHx
    Ksr = Fz * ( PKX1+  PKX2 * dfz) * np.exp( PKX3 * dfz)   
    C =  PCX1 
    mu = ( PDX1 +  PDX2 * dfz) * (1 -  PDX3 * gamma**2) * LMUX#gamma=slip ratio
    D = mu*Fz
    B = Ksr / (C * D + 0.001)  #e=0.001
    E = ( PEX1 +  PEX2 * dfz ) * (1 - ( PEX4 * np.sign(kappax))) * LEX #+ params[6] ** dfz therm fehlt
    Sv = Fz * ( PVX1 +  PVX2 * dfz) 
   
    
    FX= D * np.sin(C * np.arctan(B * kappax - E * (B*kappax - np.arctan(B * kappax)))) 
    
    return FX,D,B,C,Ksr,E


#Fy Lateral      
def Pacejka52_lat(alpha,Fz,Fz0,PVY1, PVY2, PVY3, PVY4, PKY1, PKY2, PKY3, PKY4, PKY5, PKY6, PKY7, PHY1, PHY2, PCY1, PDY1, PDY2, PDY3, PEY1, PEY2, PEY3, PEY4, PEY5, LMUY, LEY): 
    #IA=inclination angle,  Fz0 = highest load in data, 
    IA = 0

    #converts to rad
    gamma = IA*np.pi/180
    
    dfz = ((Fz-Fz0)/Fz0)  


    Svyg = Fz * ( PVY3 +  PVY4 * dfz) * gamma
    Kyg0 = Fz * ( PKY6 +  PKY7 * dfz)
    Kya =  PKY1 * Fz0 * np.sin( PKY4 * np.arctan(Fz/(( PKY2 +  PKY5 * gamma**2) * Fz0)))/(1 +  PKY3 * gamma**2)   #skalingFaktor Fz = 1, Fz0 = F´z0 da scaling factor LFZ0 = 1
    Shy = ( PHY1 +  PHY2 * dfz) + (Kyg0 * gamma - Svyg) / (Kya + 0.001)
    alphay = alpha + Shy
    
    C =  PCY1
    mu = ( PDY1 +  PDY2 * dfz) / (1 +  PDY3 * gamma**2) *LMUY
    D = mu*Fz
    B = Kya / (C * D + 0.001)  #e=0.001

    Svy = Fz * ( PVY1 *  PVY2 * dfz) + Svyg

    j=0
    Fy = np.zeros(len(alphay))
    for j in range(0,len(alphay)):
        E = ( PEY1 +  PEY2 * dfz ) * (1 +  PEY5*gamma**2 - ( PEY3 +  PEY4 * gamma) * np.sign(alphay[j])) * LEY
    
        Fy[j] = - D * np.sin(C * np.arctan(B * alphay[j] - E * (B*alphay[j] - np.arctan(B * alphay[j])))) + Svy
    
    return Fy,D,B,C,E,Kya,mu

#Mz Aligning Torque
def Pacejka52_aligning(Fz,Fz0,alpha,Fy0):


    gamma = IA*np.pi/180
    dfz = ((Fz-Fz0)/Fz0) 
    #csa = 0.9 #Vcx / (Vc + 0.1) Formel aus dem Pacejka Buch Geschwindigkeit in X Richtung / Geschwindigkeit in Contact center

    SHy = (pa.PHY1 + pa.PHY2 * dfz)
    Kya = pa.PKY1 * Fz * np.sin(pa.PKY4 * np.arctan(Fz/((pa.PKY2) * Fz)))/(1 + pa.PKY3 * gamma**2)
    Cy = pa.PCY1
    mu = (pa.PDY1 + pa.PDY2 * dfz) / (1 + pa.PDY3 * gamma**2)
    Dy = mu*Fz
    By = Kya / (Cy * Dy + 0.001)
    SVy = Fz * (pa.PVY1 * pa.PVY2 * dfz)
    Shy = (pa.PHY1 + pa.PHY2 * dfz)
    alphay = alpha + Shy
    Ey = (pa.PEY1 + pa.PEY2 * dfz ) * (1 - (pa.PEY3 + pa.PEY3 * gamma) * np.sign(alphay)) 
    Fy= Dy * np.sin(Cy * np.arctan(By * alphay - Ey * (By*alphay - np.arctan(By * alphay)))) + SVy


    Cr = 1
    Br = (pa.QBZ9 + pa.QBZ10 * By * Cy)
    Shf = SHy + SVy / (Kya + 0.001) 
    alphar = alpha + Shf
    Dr = Fz * pa.UNLOADED_RADIUS * ((pa.QDZ6 + pa.QDZ7 * dfz) + (pa.QDZ8 + pa.QDZ9 *dfz) + (pa.QDZ10 + pa.QDZ11 * dfz) * gamma *abs(gamma)) #* csa
    Mzr0 = Dr * np.cos(Cr * np.arctan(Br * alphar))
    Ct = pa.QCZ1
    Bt = (pa.QBZ1 + pa.QBZ2 *dfz + pa.QBZ3 *dfz**2) * (1 + pa.QBZ5 * abs(gamma) + pa.QBZ6 * gamma**2)
    SHt = pa.QHZ1 + pa.QHZ2 * dfz + (pa.QHZ3 + pa.QHZ4 * dfz) * gamma
    alphat = alpha + SHt
    Et = (pa.QEZ1 + pa.QEZ2 * dfz + pa.QEZ3 * dfz**2)*(1 + (pa.QEZ4 + pa.QEZ5 * gamma)*(2 / np.pi) * np.arctan(Bt * Ct * alphat))
    Dt0 = Fz * (pa.UNLOADED_RADIUS / Fz0) * (pa.QDZ1 + pa.QDZ2 *dfz)
    Dt = Dt0 * (1 + pa.QDZ3 * abs(gamma) + pa.QDZ4 * gamma**2)
    
    t0 = Dt * np.cos(Ct * np.arctan(Bt * alphat - Et * (Bt * alphat - np.arctan(Bt*alphat)))) #*csa
    Mz0 = -t0 * Fy0 + Mzr0
    
    return Mz0


#Longitudinal Fx combined sip
def Pacejka52_long_comb(Fz,Fz0,alpha,SL,Fx0): #Fx0 = Fx aus Pacejka_long

    gamma = IA*np.pi/180
    dfz = ((Fz-Fz0)/Fz0) 

    Shxa = pa.RHX1
    Exa = pa.REX1 + pa.REX2 *dfz
    Cxa = pa.RCX1
    Bxa = (pa.RBX1 + pa.RBX3 * gamma**2) * np.cos (np.arctan(pa.RBX2 * SL))
    alphas = alpha + Shxa
    Gxa0 = np.cos(Cxa * np.arctan(Bxa * Shxa - Exa *(Bxa * Shxa -np.arctan(Bxa * Shxa))))
    Gxa = np.cos(Cxa * np.arctan(Bxa * alphas - Exa *(Bxa * alphas -np.arctan(Bxa * alphas)))) / Gxa0
    Fx = Gxa * Fx0
    
    return Fx



#Lateral Fy combined slip
def Pacejka52_lat_comb(Fz,Fz0,alpha,SL,Fy0): #Fy0 = Fy aus Pacejka_lat

    gamma = IA*np.pi/180
    dfz = ((Fz-Fz0)/Fz0)

    muy = (pa.PDY1 + pa.PDY2 * dfz) / (1 + pa.PDY3 * gamma**2)
    DVyk = muy * Fz * (pa.RVY1 + pa.RVY2 * dfz + pa.RVY3 * gamma) * np.cos(np.arctan(pa.RVY4 * alpha))
    SVyk = DVyk * np.sin(pa.RVY5 * np.arctan(pa.RVY6 * SL))
    Shyk = pa.RHY1 * pa.RHY2 * dfz
    Eyk = pa.REY1 + pa.REY2 * dfz
    Cyk = pa.RCY1
    Byk = (pa.RBY1 + pa.RBY4 * gamma**2) * np.cos(np.arctan(pa.RBY2 * (alpha - pa.RBY3)))
    kappas = SL + Shyk
    Gyk0 = np.cos(Cyk * np.arctan(Byk * Shyk - Eyk *(Byk * Shyk - np.arctan(Byk * Shyk))))
    Gyk = np.cos(Cyk * np.arctan(Byk * kappas - Eyk *(Byk * kappas - np.arctan(Byk * kappas)))) / Gyk0
    Fy = Gyk * Fy0 + SVyk

    return Fy

#Alinging Torque combined slip
def Pacejka52_alinging_comb(Fz,Fz0,alpha,SL,Fx,Fy):

    gamma = IA*np.pi/180
    dfz = ((Fz-Fz0)/Fz0)
    csa = 0.9 #Vcx / (Vc + 0.1) Formel aus dem Pacejka Buch Geschwindigkeit in X Richtung / Geschwindigkeit in Contact center

    SHt = pa.QHZ1 + pa.QHZ2 * dfz + (pa.QHZ3 + pa.QHZ4 * dfz) * gamma
    alphat = alpha + SHt
    Kya = pa.PKY1 * Fz * np.sin(pa.PKY4 * np.arctan(Fz/((pa.PKY2) * Fz)))/(1 + pa.PKY3 * gamma**2)
    SVy = Fz * (pa.PVY1 * pa.PVY2 * dfz)
    Shy = (pa.PHY1 + pa.PHY2 * dfz)
    Shf = Shy + SVy / (Kya + 0.001) 
    alphar = alpha + Shf

    Kxk = Fz * (pa.PKX1+ pa.PKX2 * dfz) * np.exp(pa.PKX3 * dfz)   
    Kya = pa.PKY1 * Fz * np.sin(pa.PKY4 * np.arctan(Fz/((pa.PKY2) * Fz)))/(1 + pa.PKY3 * gamma**2)
    alphareq = np.sqrt(alphar**2 + (Kxk / Kya)**2 * SL**2) *np.sign(alphar)
    alphateq = np.sqrt(alphat**2 + (Kxk / Kya)**2 * SL**2) *np.sign(alphat)
    s = pa.UNLOADED_RADIUS * (pa.SSZ1 + pa.SSZ2 (Fy / Fz0) + (pa.SSZ3 + pa.SSZ4 * dfz)*gamma)
    Cy = pa.PCY1
    Dy = muy*Fz
    By = Kya / (Cy * Dy + 0.001) 
    muy = (pa.PDY1 + pa.PDY2 * dfz) / (1 + pa.PDY3 * gamma**2)
    DVyk = muy * Fz * (pa.RVY1 + pa.RVY2 * dfz + pa.RVY3 * gamma) * np.cos(np.arctan(pa.RVY4 * alpha))
    SVyk = DVyk * np.sin(pa.RVY5 * np.arctan(pa.RVY6 * SL))
    Dr = Fz * pa.UNLOADED_RADIUS * ((pa.QDZ6 + pa.QDZ7 * dfz) + (pa.QDZ8 + pa.QDZ9 *dfz) + (pa.QDZ10 + pa.QDZ11 * dfz) * gamma *abs(gamma)) * csa
    Cr = 1
    Br = (pa.QBZ9 + pa.QBZ10 * By * Cy)
    Dt0 = Fz * (pa.UNLOADED_RADIUS / Fz0) * (pa.QDZ1 + pa.QDZ2 *dfz)
    Dt = Dt0 * (1 + pa.QDZ3 * abs(gamma) + pa.QDZ4 * gamma**2)
    Ct = pa.QCZ1
    Bt = (pa.QBZ1 + pa.QBZ2 *dfz + pa.QBZ3 *dfz**2) * (1 + pa.QBZ5 * abs(gamma) + pa.QBZ6 * gamma**2)
    Et = (pa.QEZ1 + pa.QEZ2 * dfz + pa.QEZ3 * dfz**2)*(1 + (pa.QEZ4 + pa.QEZ5 * gamma)*(2 / np.pi) * np.arctan(Bt * Ct * alphat))
    Mzr = Dr * np.cos(Cr * np.arctan(Br * alphareq))
    Fys = Fy - SVyk
    t = Dt * np.cos(Ct * np.arctan(Bt * alphateq - Et (Bt * alphateq - np.arctan(Bt * alphateq)))) * csa
    Mzs = -t *Fys
    Mz = Mzs + Mzr + s *Fx

    return Mz


#Pacejka 5 Parameter

#Fx

def Pacejka5_long(SL):

    Fx= pa.PacLong_D * np.sin(pa.PacLong_C * np.arctan(pa.PacLong_D * SL - pa.PacLong_E * (pa.PacLong_B*SL - np.arctan(pa.PacLong_B * SL)))) 
    
    return Fx
    
#Fy
def Pacejka5_lat(alpha):

    Fy= pa.PacLat_D * np.sin(pa.PacLat_C * np.arctan(pa.PacLat_D * alpha - pa.PacLat_E * (pa.PacLat_B*alpha - np.arctan(pa.PacLat_B * alpha)))) 
    
    return Fy
    
  
    
    
    
    
    
    
    
    