UNLOADED_RADIUS     = 0.2059        # [  m ] Unloaded radius
        
## [VERTICAL]
FZ1                  = 650        # [  N ] Reference wheel load 1
FZ2                  = 1200      # [  N ] Reference wheel load 2
CZ1                  = 88029        # [ N/m] Tire vertical stiffness @FZ1
CZ2                  = 96545        # [ N/m] Tire vertical stiffness @FZ2
DZ                   = 0        # [Ns/m] Tire vertical damping
LR1                  = 0.3500        # [  - ] Coefficient for effective rolling radius @FZ1
LR2                  = 0.3500        # [  - ] Coefficient for effective rolling radius @FZ2
        
## [VERTICAL_FORCE_RANGE]
FZMIN                = 0        # [  N ] Minimum allowed wheel load
FZMAX                = 100000 #inf      # [  N ] Maximum allowed wheel load
        
## [LONG_SLIP_RANGE]
KPUMIN               = -1.5     # [  - ] Minimum valid wheel slip
KPUMAX               =  1.5     # [  - ] Maximum valid wheel slip
        
     ## [SLIP_ANGLE_RANGE]
ALPMIN               = -1.5708  # [ rad] Minimum valid slip angle
ALPMAX               =  1.5708  # [ rad] Maximum valid slip angle
        
## [INCLINATION_ANGLE_RANGE]
CAMMIN               = -0.20996 # [ rad] Minimum valid camber angle
CAMMAX               =  0.20996 # [ rad] Maximum valid camber angle
        
 ## [LONGITUDINAL_COEFFICIENTS]
DFX1                 = 32529       # [ N/#] Longitudinal slip stiffness @FZ1
DFX2                 = 51543                 # [ N/#] Longitudinal slip stiffness @FZ2
FMX1                 = 1851.8        # [  N ] Maximum longitudinal force @FZ1
FMX2                 = 3085.8              # [  N ] Maximum longitudinal force @FZ2
SMX1                 = 0.2183                # [  - ] Longitudinal slip for FMX1
SMX2                 = 0.2324              # [  - ] Longitudinal slip for FMX2
FSX1                 = 1188.4        # [  N ] Longitudinal force for pure sliding @FZ1
FSX2                 = 1798.6                  # [  N ] Longitudinal force for sliding @FZ2
SSX1                 = 0.7725                # [  - ] Longitudinal slip for pure sliding @FZ1
SSX2                 = 0.9554             # [  - ] Longitudinal slip for pure sliding @FZ2
BRAS                 = 1        # [  - ] Factor for braking asymmetry
        
## [LATERAL_COEFFICIENTS]
DFY1                 = 30163        # [ N/#] Lateral slip stiffness @FZ1
DFY2                 = 36879       # [ N/#] Lateral slip stiffness @FZ2
FMY1                 = 1604.2       # [  N ] Maximum lateral force @FZ1
FMY2                 = 2790.1       # [  N ] Maximum lateral force @FZ2
SMY1                 = 0.2031        # [  - ] Lateral slip for FMX1
SMY2                 = 0.2305        # [  - ] Lateral slip for FMX2
FSY1                 = 1072.7        # [  N ] Lateral force for pure sliding @FZ1
FSY2                 = 1478.7        # [  N ] Lateral force for sliding @FZ2
SSY1                 = 0.7912        # [  - ] Lateral slip for pure sliding @FZ1
SSY2                 = 0.7978       # [  - ] Lateral slip for pure sliding @FZ2
            
## [ROLLING_COEFFICIENTS]
KR0                  = 0.0527        # [  - ] Rolling resistance torque coefficient
KR1                  = 0.0013        # [  - ] Rolling resistance torque depending on speed