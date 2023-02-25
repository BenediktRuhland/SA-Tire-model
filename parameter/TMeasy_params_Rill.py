#tire properties
r0 = 0.293 #unleaded radius
fz0 = 3500 #payload
cz = 190000 #vertical stiffnes

#Parameters
dfx0 = 32529
dfy0 = 30163   # init slopes in N/-
fxm  = 1851.8
fym  = 1604.2   # maximum forces in N
sxm  = 0.2183
sym  = 0.2031   # sm where f(sm) = fm
fxs  = 1188.4
fys  = 1072.7   # sliding forces in N
sxs  = 0.7725
sys  = 0.7978   # ss where f(ss) = fs 

#dynamic tire offset parameters
n0  = 0.180  # normalized tire offset n/L @ sy = 0
sy0 = 0.190  # sy where n/L passes sy-axis
syE = 0.350  # sy where n/L approaches sy-axis again
