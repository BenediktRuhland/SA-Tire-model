####
#Brush Tire Model
####

import numpy as np

#----------------------------------------------------------------
#lonitudinal Fx
#----------------------------------------------------------------

#Parameters
# Fz = normal Force [N]
# mu = peak road tire fricton coefficient
# lt = contact patch length [m]
# kt = tire tangential stiffness [N/m^2]
# lam = intial deformation

def brush_long(Sx,Fz,mu,lt,kt,lam):

    #critical slip
    Sx_crit = mu * Fz / (kt * lt * (lt + lam))

    Fx = np.zeros(len(Sx))

    for i in range(0,len(Sx)):
        if Sx[i] < Sx_crit:
            Fx[i] = kt * lt * (lam + lt / 2) * Sx[i]
            lc[i] = lt
        else:
            lc[i] = mu * Fz / (lt * kt * Sx[i]) - lam
            Fx[i] = kt * Sx[i] * (lc[i]) * (lam + lc[i] / 2) + mu * Fz * (a - lc[i] / lt)

    return Fx

#----------------------------------------------------------------
#lateral Fy
#----------------------------------------------------------------

#Parameters
# Fz = normal Force [N]
# mu = peak road tire fricton coefficient
# Ca = lateral stiffness
# lt # contact patch length

# Sy in degrees

def Brush_lat(Sy, Fz, mu, Ca, lt):

    #critical slip
    Sy_crit = mu * Fz / (2*Ca)

    Fy = np.zeros(len(Sy))
    for i in range(0, len(Sy)):
        if Sy[i] < Sy_crit:
            Fy[i] = Ca * Sy[i]
        else:
            Fy[i] = mu * Fz - (mu * Fz)**2 / (4 * Ca * np.tan(Sy[i]))
    
    return Fy
