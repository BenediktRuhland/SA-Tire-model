####
#Brush Tire Model
####

import numpy as np

import parameter.Pacejka_Params_Indy as pa
#----------------------------------------------------------------
#lonitudinal Fx
#https://skill-lync.com/student-projects/Longitudinal-brush-tire-model-50471
#----------------------------------------------------------------

#Parameters
# Fz = normal Force [N]
# mu = peak road tire fricton coefficient
# lt = contact patch length [m]
# kt = tire tangential stiffness [N/m^2]
# lam = intial deformation

def brush_long(Sx,Fz,mu,kt,lam):

    #contact patch length
    lt = np.sqrt((pa.UNLOADED_RADIUS**2) - (pa.TYRE_RADIUS_MOD**2))

    Sx = Sx * 10 #scaling factor
    #critical slip
    Sx_crit = mu * Fz / (kt * lt * (lt + lam))

    Fx = np.zeros(len(Sx))
    lc = np.zeros(len(Sx)) 

    
    for i in range(0,len(Sx)):
        
        if Sx[i] >= 0:
            if Sx[i] < Sx_crit:
                Fx[i] = kt * lt * (lam + lt / 2) * Sx[i]
                lc[i] = lt
            else:
                if (lt * kt * Sx[i]) - lam == 0:
                    lc[i] = 0
                else:
                    lc[i] = mu * Fz / (lt * kt * Sx[i]) - lam
                    Fx[i] = kt * Sx[i] * (lc[i]) * (lam + lc[i] / 2) + mu * Fz * (1 - lc[i] / lt)
        else:
            if abs(Sx[i]) < Sx_crit:
                Fx[i] = - kt * lt * (lam + lt / 2) * abs(Sx[i])
                lc[i] = lt
            else:
                lc[i] =  mu * Fz / (lt * kt * abs(Sx[i])) - lam
                Fx[i] = - (kt * abs(Sx[i]) * (lc[i]) * (lam + lc[i] / 2) + mu * Fz * (1 - lc[i] / lt))
    return Fx

#----------------------------------------------------------------
#lateral Fy
#----------------------------------------------------------------

#Parameters
# Fz = normal Force [N]
# mu = peak road tire fricton coefficient
# Ca = lateral stiffness [N /rad]
# lt # contact patch length

# Sy in degrees

def brush_lat(Sy, Fz, mu, kt, lam):

    #critical slip
    #Sy_crit = mu * Fz / (2*Ca)

    #Fy = np.zeros(len(Sy))
    #for i in range(0, len(Sy)):

    #    if Sy[i] >= 0:
    #        if Sy[i] < Sy_crit:
    #            Fy[i] = Ca * Sy[i]
    #        else:
    #            Fy[i] = mu * Fz - (mu * Fz)**2 / (4 * Ca * np.tan(Sy[i]))
    #    else:
    #        if abs(Sy[i]) < Sy_crit:
    #            Fy[i] = - Ca * abs(Sy[i])
    #        else:
    #            Fy[i] = - (mu * Fz - (mu * Fz)**2 / (4 * Ca * np.tan(abs(Sy[i]))))
    #return Fy

    #contact patch length
    lt = np.sqrt((pa.UNLOADED_RADIUS**2) - (pa.TYRE_RADIUS_MOD**2))
 
    Sy = Sy * 10 #scaling factor
    #critical slip
    Sx_crit = mu * Fz / (kt * lt * (lt + lam))

    Fx = np.zeros(len(Sy))
    lc = np.zeros(len(Sy)) 

    
    for i in range(0,len(Sy)):
        
        if Sy[i] >= 0:
            if Sy[i] < Sx_crit:
                Fx[i] = kt * lt * (lam + lt / 2) * Sy[i]
                lc[i] = lt
            else:
                if (lt * kt * Sy[i]) - lam == 0:
                    lc[i] = 0
                else:
                    lc[i] = mu * Fz / (lt * kt * Sy[i]) - lam
                    Fx[i] = kt * Sy[i] * (lc[i]) * (lam + lc[i] / 2) + mu * Fz * (1 - lc[i] / lt)
        else:
            if abs(Sy[i]) < Sx_crit:
                Fx[i] = - kt * lt * (lam + lt / 2) * abs(Sy[i])
                lc[i] = lt
            else:
                lc[i] =  mu * Fz / (lt * kt * abs(Sy[i])) - lam
                Fx[i] = - (kt * abs(Sy[i]) * (lc[i]) * (lam + lc[i] / 2) + mu * Fz * (1 - lc[i] / lt))
    return Fx

print(np.sqrt((pa.UNLOADED_RADIUS**2) - (pa.TYRE_RADIUS_MOD**2)))