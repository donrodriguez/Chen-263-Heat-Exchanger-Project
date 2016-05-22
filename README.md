# Chen-263-Heat-Exchanger-Project

from numpy import *
import scipy as spy


print("Substances available for this heat exchanger calculation: ")
subs = ["Water", "R134a", "Ethanol", "2,2,4-trimethylpentane (tmp for short)"]
for x in subs:
    print(x)

#Unit conversions from AES to SI

def U_SI(x): #Converts Overall Heat Transfer Coefficient (U) in AES units (Btu/(hr-ft^2-F)) to SI
    return x* 3.2808**2 * 1.8 /(3600 * 9.486*10**(-4))
def mflowSI(x): #converts lb-m/s to kg/s
    return x*.453593
def tempSI(F): #Converts F to C degrees
    return (F-32)/1.8

water = {"MW": 18.01528, "Mp": 273.15, "Bp" : 373.15 , "A": 276370 , "B": -2090.1 , "C":8.125, "D":-0.014116,
         "E":9.3701E-06}
r134a = {"MW":102.03089, "Mp":172.00 , "Bp" :247.08 , "A":6.5108E+05 , "B":-9.5057E+03 , "C":6.2835E+01 ,
         "D":-1.8264E-01, "E":2.0031E-04 }
ethanol = {"MW":46.06844 , "Mp":159.05, "Bp":351.44 , "A":1.0264E+05, "B":-1.3963E+02, "C":-3.0341E-02,
           "D":2.0386E-03 , "E":0.0 }
tmp = {"MW":114.22852 , "Mp":165.777, "Bp":372.388, "A":9.5275E+04, "B":6.9670E+02 , "C":-1.3765E+00, "D":2.1734E-03,
       "E":0.0}

Unit_type = input("What units will you be using (AES or SI) for your input calculations?:  ")
if Unit_type.upper() == "AES":
    hotf = input("Enter the hotter fluid species: ")



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
