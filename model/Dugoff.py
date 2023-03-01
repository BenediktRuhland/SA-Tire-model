#----------------------------------------------------------------
#Dugoff Tire Model

#Dugoff Modell has no peak force --> modified Dugoff model (skaling factir G added)
#Paper: A Dynaic Model for Tire/road Friction under Combined Longitudina/Lateral Slip Situation, SAE International, 2014
#----------------------------------------------------------------

import numpy as np

#Longitudenal Fx

def dugoff_long(Sx,Sy,Fz,Cs,Ca,mu):

    Fx = np.zeros(len(Sx))
    for i in range(0,len(Sx)):

        if (2 * np.sqrt((Cs * Sx[i])**2 + (Ca * np.tan(Sy[i]))**2 )) == 0:
            lam = 0
        else:
            lam = (mu * Fz * (1 - Sx[i])) / (2 * np.sqrt((Cs * Sx[i])**2 + (Ca * np.tan(Sy[i]))**2 ))

        if lam < 1:
            f_lam = (2 - lam) * lam
        else:
            f_lam = 1

        Gx =  (1.15 - 0.75 * mu)  * abs(Sx[i])**2 - (1.63 - 0.75 * mu) * abs(Sx[i]) + 1.5

        if (1 - Sx[i]) == 0:
            Fx[i] = Fx[i-1]
        else:
            Fx[i] = Cs * Sx[i] * f_lam * Gx / (1 - Sx[i])

    return Fx

#lateral Fy

def dugoff_lat(Sx,Sy,Fz,Cs,Ca,mu):

    Fy = np.zeros(len(Sy))
    for i in range(0,len(Sy)):

        if (2 * np.sqrt((Cs * Sx[i])**2 + (Ca * np.tan(Sy[i]))**2 )) == 0:
            lam = 0
        else:
            lam = (mu * Fz * (1 - Sx[i])) / (2 * np.sqrt((Cs * Sx[i])**2 + (Ca * np.tan(Sy[i]))**2 ))

        if lam < 1:
            f_lam = (2 - lam) * lam
        else:
            f_lam = 1

        Gy = ((mu - 1.6) * np.tan(abs(Sy[i])) + 1.5) * 10
        if (1 - Sx[i]) == 0:
            Fy[i] = Fy[i-1]
        else:
            Fy[i] = Ca * np.tan(Sy[i]) * f_lam * Gy/ (1 - Sx[i])

    return Fy