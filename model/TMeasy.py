#######################################
#TMeasy Tire Model
#Prof.Rill Vehical Dynamics (Vorlesungsaufzeichnungen) Kapitel 3.
#TUfast TMeasy Tire Model für Lapsim 
########################################

import numpy as np

def calcFY(SX, SY, DFY, FMY, FSY, SMY, SSY,SMX):
    DFY_arr = np.ones(len(SX)) + DFY
    FMY_arr = np.ones(len(SX)) + FMY
    FSY_arr = np.ones(len(SX)) + FSY
    SMY_arr = np.ones(len(SX)) + SMY
    SSY_arr = np.ones(len(SX)) + SSY
    SMX_arr = np.ones(len(SX)) + SMX

    #Normalized slip
    SIGX = SX/SMX_arr
    SIGY = SY/SMY_arr

    #Combined Slip
    SIG = np.sqrt(SIGX**2 + SIGY**2)

    Fy = np.zeros(len(SIG)) #np.zeros(a)
    j=0
    for j in range(0,len(SIG)):
        if SIG[j] <= 1:
            AY1 = DFY_arr[j]*SMY_arr[j]
            BY1 = SMY_arr[j]/FMY_arr[j]*DFY_arr[j]-2
            Fy[j] = AY1*SIGY[j]/(SIG[j]**2+BY1*SIG[j]+1)
        elif SIG[j] <= SSY_arr[j]/SMY_arr[j]:
            AY2 = (SIG[j]-1)/(SSY_arr[j]/SMY_arr[j]-1)
            Fy[j] = SIGY[j]/SIG[j] *(FMY_arr[j] - (FMY_arr[j]-FSY_arr[j])*AY2**2*(3-2*AY2))
        else:
            Fy[j] = SIGY[j]/SIG[j] * FSY_arr[j]
    return Fy

def calcFx( SX, SY, DFX, FMX, FSX, SMX,SSX, SMY):
    DFX_arr = np.ones(len(SX)) * DFX
    FMX_arr = np.ones(len(SX)) * FMX
    FSX_arr = np.ones(len(SX)) * FSX
    SMX_arr = np.ones(len(SX)) * SMX
    SSX_arr = np.ones(len(SX)) * SSX
    SMY_arr = np.ones(len(SX)) * SMY

    #Normalized slip
    SIGX = SX/SMX_arr
    SIGY = SY/SMY_arr

    #Combined Slip
    SIG = np.sqrt(SIGX**2 + SIGY**2)

    FX = np.zeros(len(SIG))
    j=0
    for j in range(0,len(SIG)):
        if SIG[j]<=1:
            AX1 = DFX_arr[j]*SMX_arr[j]
            BX1 = SMX_arr[j]/FMX_arr[j]*DFX_arr[j]-2
            FX[j] = AX1*SIGX[j]/(SIG[j]**2+BX1*SIG[j]+1)
        elif SIG[j] <= SSX_arr[j]/SMX_arr[j]:
            AX2 = (SIG[j]-1)/(SSX_arr[j]/SMX_arr[j]-1)
            FX[j] = SIGX[j]/SIG[j] *(FMX_arr[j] - (FMX_arr[j]-FSX_arr[j])*AX2**2*(3-2*AX2))
        else:
            FX[j] = SIGX[j]/SIG[j] * FSX_arr[j]
            
    return FX


