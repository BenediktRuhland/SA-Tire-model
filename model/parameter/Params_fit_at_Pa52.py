#----------------------------------
#Parameter verschiedener Reifenmodelle gefitted auf Pacejka52 mit Pacejka_params_Indy
#Pacejka52 orginal Kurve
#----------------------------------

#TMeasy

#Fx: 
DFX_long = 2.73698918e+05
FMX_long = 6.49227443e+03
FSX_long = 4.92427382e+03
SMX_long = 5.65704982e-02
SSX_long = 2.59229223e-01
SMY_long = 5.03325089e+00

#Fy: 
#DFY_lat = 9.01646321e+04
#FMY_lat = 6.39266085e+03
#FSY_lat = 5.30979506e+03
#SMY_lat = 1.15048404e-01
#SSY_lat = 3.01495645e-01
#SMX_lat = 3.54195854e+00

#Fy neu:
DFY_lat = 9.46557086e+04  
FMY_lat = 6.40338330e+03  
FSY_lat = 5.57015421e+03 
SMY_lat = -8.82115882e-01
SSY_lat = -7.30516695e-01  
SMX_lat = 2.03471709e+00
#MFsimple

#Fx
B_long =19.7526835
C_long =1.99999146
D_long =1.63455163
E_long = 0.75471735

#Fy
B_lat = 9.10763555
C_lat = 1.5991035
D_lat = 1.60244708 
E_lat = -1.89562919
   

#Brush

#Fx
mu_b_long  = 1.62344260e+00
kt_b_long  = 3.63683914e+06
lam_b_long = 3.10083340e-02

#Fy
mu_b_lat = 1.58991258e+00
Ca_b_lat = 1.05121281e+06
lam_b_lat = 5.74790914e-02

#Dugoff

#Fx
Cs_long = 1.58443682e+05
Ca_long = -1.59200977e+00
mu_long = 1.24651131e+00
#mit Syyy = np.linspace(0.001,0.2,1000)

#Fy
Cs_lat = 2.74289218e-01
Ca_lat = 6.16477615e+03
mu_lat =  1.43548745e-01
#mit Sxxx = np.linspace(0,0.1,1000)
#Fiala

#Fy
fi_C_lat = -1.26902260e+05
fi_mu1 = -2.32752807e+00
fi_mu2 = -1.44561563e+00
