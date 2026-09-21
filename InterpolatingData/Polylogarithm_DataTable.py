#
#                                        Polylogarithm_DataTable.py
#                                           31 Jul. 2025
#                           Generates data table of polylogarithm functions 
#
#   Li_{-1/2}(-exp(z)), Li_{1/2}(-exp(z)), Li_{3/2}(-exp(z)), Li_{5/2}(-exp(z)), Li_{7/2}(-exp(z))
#
#       as a function of z. Outputs data table "Polylogarithms.dat" used in "InterpolatorGeneration.cpp"
#       which generates interpolating functions used in "ElectronMHDAxisymmetricSpherical" C++ project
#
#########################################################################################################

import numpy as np
from mpmath import mp

PolyLogN1_2Vec = np.frompyfunc(lambda *a: float( mp.re(mp.polylog(-1/2,*a)) ), 1, 1)
PolyLog1_2Vec = np.frompyfunc(lambda *a: float( mp.re(mp.polylog(1/2,*a)) ), 1, 1)
PolyLog3_2Vec = np.frompyfunc(lambda *a: float( mp.re(mp.polylog(3/2,*a)) ), 1, 1)
PolyLog5_2Vec = np.frompyfunc(lambda *a: float( mp.re(mp.polylog(5/2,*a)) ), 1, 1)
PolyLog7_2Vec = np.frompyfunc(lambda *a: float( mp.re(mp.polylog(7/2,*a)) ), 1, 1)
PolyLog9_2Vec = np.frompyfunc(lambda *a: float( mp.re(mp.polylog(9/2,*a)) ), 1, 1)

num = 400 #number of data table entries

z_ar = np.linspace(-50.5,50.5,num,endpoint=True)

PLN1_2 = PolyLogN1_2Vec(-np.exp(z_ar))
PL1_2 = PolyLog1_2Vec(-np.exp(z_ar))
PL3_2 = PolyLog3_2Vec(-np.exp(z_ar))
PL5_2 = PolyLog5_2Vec(-np.exp(z_ar))
PL7_2 = PolyLog7_2Vec(-np.exp(z_ar))

# z_ar = np.logspace(-22,22,num,endpoint=True)

# PLN1_2 = PolyLogN1_2Vec(-z_ar)
# PL1_2 = PolyLog1_2Vec(-z_ar)
# PL3_2 = PolyLog3_2Vec(-z_ar)
# PL5_2 = PolyLog5_2Vec(-z_ar)
# PL7_2 = PolyLog7_2Vec(-z_ar)

with open("Polylogarithms.dat", 'w') as file:
    file.write("z   Li_{-1/2}   Li_{1/2}   Li_{3/2}   Li_{5/2}   Li_{7/2}\n")
    for z, PLN1_2, PL1_2, PL3_2, PL5_2, PL7_2 in zip(z_ar, PLN1_2, PL1_2, PL3_2, PL5_2, PL7_2):
        file.write(f"{z}   {PLN1_2}   {PL1_2}   {PL3_2}   {PL5_2}   {PL7_2}\n")
