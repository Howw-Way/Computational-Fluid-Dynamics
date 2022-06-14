# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 10:23:37 2020

@author: Hao
"""
#This algorithm use Davis mechanism to calculate laminar flame speed,then use limit flame method to calculate the upper limit of C3H8
import cantera as ct
import numpy as np
import pandas as pd
import time
start = time.clock()

width = 0.057  # m
loglevel = 1  # amount of diagnostic output (0 to 8)
a=1.11
O2=1
Ar=0.25
CO2=0
N2=0
p = 1*101325  # pressure [Pa]
Tin =450.0  # unburned gas temperature [K]

while 1:
# Simulation parameters
    C3H8=a
    reactants = 'C3H8:{}, O2:{}, CO2:{}, Ar:{},N2:{}' .format(C3H8,O2,CO2,Ar,N2) 


# IdealGasMix object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
    gas = ct.Solution('DavischemAr.CTI')
    gas.TPX = Tin, p, reactants
    Udensity=gas.density
    Cp=gas.cp
    y=gas.thermal_conductivity
    
# Set up flame object
    f= ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
    #f.solve(loglevel=loglevel, auto=True)
    f.show_solution()

    
# Solve with mixture-averaged transport model
    f.transport_model = 'Mix'
    f.solve(loglevel=loglevel, auto=True)

# Solve with the energy equation enabled
    f.save('c3h8_adiabatic.xml', 'mix', 'solution with mixture-averaged transport')
    f.show_solution()
    
# Solve with multi-component transport properties
    f.transport_model = 'Multi'
    f.solve(loglevel) # don't use 'auto' on subsequent solves
    f.show_solution()
    Su=f.u[0]
    f.save('c3h8_adiabatic.xml','multi', 'solution with multicomponent transport')
    
# write the velocity, temperature, density, and mole fractions to a CSV file
    f.write_csv('c3h8_adiabatic.csv', quiet=False)
    
#Open CSV fiel to get density
    with open("c3h8_adiabatic.csv") as csvfile:
        mLines = csvfile.readlines()
    targetLine = mLines[-2]
    Mdestiny=targetLine.split(',')[4]
    Tflame=targetLine.split(',')[3]
#print('Udestiny = {0:7f} kg/m3'.format(eval(Mdestiny)))
    Bdestiny=eval(Mdestiny)

#Su limt
    x=y/(Udensity*Cp)
    flim=((2*9.8*Bdestiny*y)/(Cp*Udensity))**(1/3)
#print result
    m=100*(a/(a+O2+CO2+N2+Ar))
    print('m= {0:7f} %'.format(m))
    print('a= {0:7f} '.format(a))
    print('Su= {0:7f} m/s'.format(Su))
    print('flim = {0:7f} m/s'.format(flim))
    print('Tflame= {0:7f} K'.format(eval(Tflame)))
    elapsed = (time.clock() - start)
    print("Time used:",elapsed)
#entering loop to check the condition
    if Su<0.49*flim:
        a=a-0.001
    elif Su>0.51*flim:
        a=a+0.001
    else:
        break