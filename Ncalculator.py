#!/usr/bin/env python3

# Program to calculate the eclipse number.
# It prompts for the Eclipse time and error, and the T(0) and error
# It returns the eclipse number and error.

import math

def N(T,T0,P):
    return (T-T0)/P

def N_err(T,T0,P,T_eff,T0_eff,P_eff):
    return N*math.sqrt((math.sqrt(T_err**2+T0_err**2)/(T-T0))**2+(P_err/P)**2)

T = float(input ("T: "))
T_err = float(input ("T_err: "))
T0 = float(input ("T0: "))
T0_err = float(input ("T0_err: "))
P = float(input ("P: "))
P_err = float(input ("P_err: "))
N = N(T,T0,P)
N_err = N_err(T,T0,P,T_err,T0_err,P_err)

print ("N: ", N)
print ("N_err: ", N_err)