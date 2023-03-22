####
#Fiala Tire Model
#
####y

import numpy as np

#longtudinal Fx

def fiala_long(Sx,Fz,C_long,mu1,mu2):

    #Matlab dokumentation https://de.mathworks.com/help/vdynblks/ref/fialawheel2dof.html
    
    #mu1 = Coefficient of static friction
    #mu2 = Coefficient of kinetic friction
    #Fx = np.zeros(len(Sx))
    
    #for i in range(len(Sx)):
    #    kappa = np.sqrt(Sx[i]**2 + np.tan(Sy[i]**2))
    #    lam = 1 #scaling factor frictiom
    #    mu = (mu1 - ( mu1 - mu2) * kappa) * lam 
    #    Sx_crit = abs((mu1 * Fz) / (2 * C_long))
        
    #    if abs(Sx[i]) < Sx_crit:
    #        Fx[i] = C_long * Sx[i]
    #    else:
    #        Fx[i] = np.tanh(4 * Sx[i]) * (mu * Fz - abs((mu * Fz)**2 / (4 * Sx[i] *C_long)))
        
    #return Fx ,Sx_crit

    Sy_crit = np.arctan(3 * mu1 * Fz / C_long)
    
    a1 = -C_long
    a2 = ((C_long**2) / (3 * mu1 * Fz)) * (2 - (mu2/mu1))
    a3 = -((C_long**3) / (9 * mu1**2 * Fz **2)) * (1- ((2 * mu2) / (3 * mu1)))

    Fy = np.zeros(len(Sx))

    for i in range(0, len(Sx)):
        if abs(Sx[i]) < Sy_crit:
            Fy[i] = a1 * np.tan(Sx[i]) + a2 * abs(np.tan(Sx[i])) * np.tan(Sx[i]) + a3 * np.tan(Sx[i])**3
        else:
            Fy[i] = -mu2 * Fz * np.sign(Sx[i]) 
    return Fy

def fiala_lat(Sy,Fz,C_lat,mu1,mu2):

    #Fy = np.zeros(len(Sy))
    #for i in range(len(Sy)):

    #    beta = np.sqrt(Sx[i]**2 + np.tan(Sy[i])**2)
    #    if (mu2 - beta * (mu2 - mu1)) > 0:
    #        mu = (mu2 - beta * (mu2 - mu1))
    #    else:
    #        mu = 0

    #    Sy_crit = abs(np.arctan(3 * mu * Fz / C_lat))
    #    H = 1 - (1/3) * C_lat * abs(np.tan(Sy)) / (mu * Fz)
    #    inside = (mu1 * Fz)**2 
    #    inside = max(0,inside)
    #    zeta = np.sqrt(inside) / (mu1 * Fz)
    #    if abs(Sy[i]) < Sy_crit:
    #        #Fy[i] = -mu * Fz * sign(Sy[i]) * (1 - H**3)
    #        Fy[i] = -C_lat * np.tan(Sy[i]) + (C_lat**2/(3 * zeta * mu1*Fz)) * abs(np.tan(Sy[i]))*np.tan(Sy[i]) - (C_lat**3/(27*zeta**2*mu1**2*Fz**2))*np.tan(Sy[i])**3
    #    else:
    #        Fy[i] = -mu * Fz * np.sign(Sy[i]) 

    Sy_crit = np.arctan(3 * mu1 * Fz / C_lat)
    
    a1 = -C_lat
    a2 = ((C_lat**2) / (3 * mu1 * Fz)) * (2 - (mu2/mu1))
    a3 = -((C_lat**3) / (9 * mu1**2 * Fz **2)) * (1- ((2 * mu2) / (3 * mu1)))

    Fy = np.zeros(len(Sy))

    for i in range(0, len(Sy)):
        if abs(Sy[i]) < Sy_crit:
            Fy[i] = a1 * np.tan(Sy[i]) + a2 * abs(np.tan(Sy[i])) * np.tan(Sy[i]) + a3 * np.tan(Sy[i])**3
        else:
            Fy[i] = -mu2 * Fz * np.sign(Sy[i]) 
    return Fy