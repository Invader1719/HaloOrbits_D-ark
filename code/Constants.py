from  Functions import Solver5
import numpy as np


class Const:
    Mass_Earth = 5.9722e24 #kg
    Mass_Moon = 7.3477e22 #kg
    Mass_Sun = 1.988416e30 #kg
    Grav = 6.67430e-11 #CI units
    Dist_Moon =  384399000 #m
    cos30 =( 3 ** 0.5) / 2
    cos45 =  1 / (2 ** 0.5)

    #Обезразмереннные величины

    Mu = Mass_Moon/(Mass_Earth + Mass_Moon)
    

def countlibra(GivenMu):
    L1 = np.array([1 - GivenMu - float(Solver5([1, -(3 - GivenMu), (3 - 2 * GivenMu), - GivenMu, + (2 * GivenMu), -GivenMu]).RealRoots()[0]), 0, 0])
    L2 = np.array([1 - GivenMu + float(Solver5([1, (3 - GivenMu), (3 - 2 * GivenMu), - GivenMu, -2 * GivenMu, -GivenMu]).RealRoots()[0]), 0, 0])
    L3 = np.array([-GivenMu - 7/12 * GivenMu, 0, 0])
    L4 = np.array([0.5 - GivenMu, Const.cos30, 0])
    L5 = np.array([0.5 - GivenMu, - Const.cos30, 0])

    return [[L1[0], L2[0], L3[0], L4[0], L5[0]],
                 [L1[1], L2[1], L3[1], L4[1], L5[1]],
                 [L1[2], L2[2], L3[2], L4[2], L5[2]]]
        
