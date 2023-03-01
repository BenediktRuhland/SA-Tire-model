####
#Fiala Tire Model
#
####

import numpy as np

#longtudinal Fx

def fiala_long(Sx,Sy,Fz,C_long,mu1,mu2):

    Fx = np.zeros(len(Sx))
    for i in range(len(Sx)):

        beta = np.sqrt(Sx[i]**2 + np.tan(Sy[i])**2)
        if (mu2 - beta * (mu2 - mu1)) > 0:
            mu = (mu2 - beta * (mu2 - mu1))
        else:
            mu = 0

        Sx_crit = 0.5 * abs((mu * Fz) / C_long)

        if abs(Sx[i]) < Sx_crit:
            Fx[i] = C_long * Sx[i]
        else:
            Fx[i] = (mu * Fz - (mu * Fz / 2)**2) / C_long * Sx[i]**2
    
    return Fx


#lateral Fy

def fiala_lat(Sx,Sy,Fz,C_lat, mu1, mu2):

    Fy = np.zeros(len(Sy))
    for i in range(len(Sy)):

        beta = np.sqrt(Sx[i]**2 + np.tan(Sy[i])**2)
        if (mu2 - beta * (mu2 - mu1)) > 0:
            mu = (mu2 - beta * (mu2 - mu1))
        else:
            mu = 0

        Sy_crit = np.arctan(3 * mu * Fz / C_lat)
        H = 1 - (1/3) * C_lat * abs(np.tan(Sy)) / (mu * Fz)

        if abs(Sy[i]) < Sy_crit:
            Fx[i] = -mu * Fz * sign(Sy[i]) * (1 - H**3)
        else:
            Fx[i] = -mu * Fz * sign(Sy[i]) 

    return Fx