# Chen-263-Heat-Exchanger-Project
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
from scipy.optimize import fsolve


print("Substances available for this heat exchanger calculation: ")
subs = ["Water", "R134a", "Ethanol", "2,2,4-trimethylpentane (type tmp for short)"]
for x in subs:
    print(x)
#how does the user imput a substance?

#Unit conversions from AES to SI
#User will need to be able to imput numbers that can then be changed to SI
def U_SI(x): #Converts Overall Heat Transfer Coefficient (U) in AES units (Btu/(hr-ft^2-F)) to SI
    return x* 3.2808**2 * 1.8 /(3600 * 9.486*10**(-4))
def mflowSI(x): #converts lb-m/s to kg/s
    return x*.453593
def tempSI(F): #Converts F to C degrees
    return (F-32)/1.8

#Substance properties
water = {"MW": 18.01528, "Mp": 273.15, "Bp" : 373.15 , "A": 276370.0 , "B": -2090.1 , "C":8.125, "D":-0.014116,
         "E":9.3701E-06}
r134a = {"MW":102.03089, "Mp":172.00 , "Bp" :247.08 , "A":6.5108E+05 , "B":-9.5057E+03 , "C":6.2835E+01 ,
         "D":-1.8264E-01, "E":2.0031E-04 }
ethanol = {"MW":46.06844 , "Mp":159.05, "Bp":351.44 , "A":1.0264E+05, "B":-1.3963E+02, "C":-3.0341E-02,
           "D":2.0386E-03 , "E":0.0 }
tmp = {"MW":114.22852 , "Mp":165.777, "Bp":372.388, "A":9.5275E+04, "B":6.9670E+02 , "C":-1.3765E+00, "D":2.1734E-03,
       "E":0.0}

#Creates polynomial with Tc_o as the variable
def Tc_o_func(T, Th_i, Tc_i, Th_o, mh, mc, hf, cf):
    '''
    :param Th_i: input Temp for hot fluid
    :param Tc_i: input Temp for cold fluid
    :param Th_o: output Temp for hot fluid
    :param mh: mass flow rate of hot fluid
    :param mc: mass flow rate of cold fluid
    :param hf: hot fluid species
    :param cf: cold fluid species
    :param T: Temperature
    :return: a polynomial with Tc_o as the variable
    '''

    Th_avg = (Th_i + Th_o)/2  # average Temp of the hot fluid
    Cp_h = hf["A"] + hf["B"] * Th_avg + hf["C"] * Th_avg ** 2 + hf["D"] * Th_avg ** 3 + \
           hf["E"] * Th_avg ** 4  # heat capacity of hot fluid
    beta = cf["MW"] * mh * Cp_h * (Th_i - Th_o) / (mc*hf["MW"])
    a = cf["A"] + (1/4)*(-cf["C"]*Tc_i**2 - cf["D"] * Tc_i**3) - 3*cf["E"]*Tc_i**4/16
    b = cf["B"]/2 + cf["C"]*Tc_i/4 - cf["E"]*Tc_i**3/8
    c = cf["C"]/4 + cf["D"]*Tc_i/4 + cf["E"]*Tc_i**2/8
    d = cf["D"]/8 + 3*cf["E"]*Tc_i/16
    e = cf["E"]/16
    alpha = -cf["A"]*Tc_i - cf["B"]*Tc_i**2/2 - cf["C"]*Tc_i**3/4 - cf["D"]*Tc_i**4/8 - cf["E"]*Tc_i**5/16
    y = a*T + b*T**2 + c*T**3 + d*T**4 + e*T**5 + alpha - beta
    return y

#creates polynomial with Th_o as the variable
def Th_o_func(T, Tc_i, Th_i, Tc_o, mh, mc, cf, hf):
    '''
    :param Tc_i: input Temp for cold fluid
    :param Th_i: input Temp for hot fluid
    :param Tc_o: output Temp for cold fluid
    :param mh: mass flow rate of hot fluid
    :param mc: mass flow rate of cold fluid
    :param cf: cold fluid species
    :param hf: hot fluid species
    :param T: temperature
    :return: a polynomial with Th_o as the variable
    '''

    Th_avg = (Tc_i + Tc_o) / 2  # average Temp of the cold fluid
    Cp_c = cf["A"] + cf["B"] * Th_avg + cf["C"] * Th_avg ** 2 + cf["D"] * Th_avg ** 3 + \
           cf["E"] * Th_avg ** 4  # heat capacity of cold fluid
    beta = hf["MW"] * mc * Cp_c * (Tc_o - Tc_i) / (mh * cf["MW"])
    a = -hf["A"] + (1 / 4) * (hf["C"] * Th_i ** 2 + hf["D"] * Th_i ** 3) + 3 * hf["E"] * Th_i ** 4 / 16
    b = -hf["B"] / 2 - hf["C"] * Th_i / 4 + hf["E"] * Th_i ** 3 / 8
    c = -hf["C"] / 4 - hf["D"] * Th_i / 4 - hf["E"] * Th_i ** 2 / 8
    d = -hf["D"] / 8 - 3 * hf["E"] * Th_i / 16
    e = -hf["E"] / 16
    alpha = hf["A"] * Th_i + hf["B"] * Th_i ** 2 / 2 + hf["C"] * Th_i ** 3 / 4 + hf["D"] * Th_i ** 4 / 8 \
            + hf["E"] * Th_i ** 5 / 16
    y = a*T + b*T**2 + c*T**3 + d*T**4 + e*T**5 + alpha - beta
    return y
    
# Asks User input
hotf = input("Enter Hotter fluid species: ")

while True:
    Thi, Unit1 = input("Enter Hot Fluid inlet Temp & Units (F or K): ").split()
    if Unit1.upper()!="F" and Unit1.upper()!="K":
        print("Invalid Unit Entries.")
        continue
    try:
        Thi = float(Thi)
    except ValueError:
        print("Invalid Temp Entry")
        continue
    if Unit1.upper() == "F":
        Thi = tempSI(Thi)
    elif Unit1.upper() == "K":
        Thi = Thi
    if Thi < substances[hotf]["Mp"]:
        print("Invalid Temp. Solid substance.")
        continue
    elif Thi > substances[hotf]["Bp"]:
        print("Invalid Temp. Gaseous substance.")
        continue
    else:
        break

coldf = input("Enter Colder fluid species: ")

while True:
    Tci, Unit2 = input("Enter Cold Fluid inlet Temp & Units (F or K): ").split()
    if Unit2.upper()!="F" and Unit2.upper()!="K":
        print("Invalid Unit Entries.")
        continue
    try:
        Tci = float(Tci)
    except ValueError:
        print("Invalid Temp Entry")
        continue
    if Unit2.upper() == "F":
        Tci = tempSI(Tci)
    elif Unit2.upper() == "K":
        Tci = Tci
    if Tci < substances[coldf]["Mp"]:
        print("Invalid Temp. Solid substance.")
        continue
    elif Tci > substances[coldf]["Bp"]:
        print("Invalid Temp. Gaseous substance.")
        continue
    else:
        break
        


# import Excel file of Data
#knowns
#Equations
T1=Th,i-Tc,o
T2=Th,o-Tc,i
Tlm=(T2-T1)/ln(T2/T1)
R=(Th,i-Th,o)/(Tc,o-Tc,i)
P=(Tc,o-Tc,i)/(Th,i-Tc,i)
F=(sqrt(R**2+1)*ln[(1-P)/(1-P*R)])/((R-1)*ln((2-P*(R+1-sqrt(R**2+1)))/(2-P*(R+1+sqrt(R**2+1))))
q=F*U*A*Tlm
q=mh*Cph*(Th,i-Th,o)
q=mc*Cpc*(Tc,o-Tc,i)
#Answers

A=F*U*Tlm/q
cost=A*1000
print(A)
print('$',cost)
print(Th_o or Tc_o)#will put in equation for this
