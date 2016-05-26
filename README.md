# Chen-263-Heat-Exchanger-Project
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

# Lists available species for this program
print("Substances available for this heat exchanger calculation: ")
subs = ["Water", "R134a", "Ethanol", "2,2,4-trimethylpentane (Enter tmp)"]
for x in subs:
    print(x)
print(" ")

# Unit conversions from AES to SI

def U_SI(x):  # Converts Overall Heat Transfer Coefficient (U) in AES units (Btu/(hr-ft^2-F)) to SI
    return x * 3.2808**2 * 1.8 / (3600 * 9.486*10**(-4))
def mflowSI(x):  # converts lb-m/s to kg/s
    return x*.453593
def tempSI(F):  # Converts F to K degrees
    return (F-32)/1.8 + 273.15

# Substance properties
water = {"MW": 18.01528, "Mp": 273.15, "Bp": 373.15, "A": 276370.0, "B": -2090.1, "C": 8.125, "D": -0.014116,
         "E": 9.3701E-06}
r134a = {"MW": 102.03089, "Mp": 172.00, "Bp": 247.08, "A": 6.5108E+05, "B": -9.5057E+03, "C": 6.2835E+01,
         "D": -1.8264E-01, "E": 2.0031E-04}
ethanol = {"MW": 46.06844, "Mp": 159.05, "Bp": 351.44, "A": 1.0264E+05, "B": -1.3963E+02, "C": -3.0341E-02,
           "D": 2.0386E-03, "E": 0.0}
tmp = {"MW": 114.22852, "Mp": 165.777, "Bp": 372.388, "A": 9.5275E+04, "B": 6.9670E+02, "C": -1.3765E+00,
       "D": 2.1734E-03, "E": 0.0}
# Dictionary of dictionaries
subst = {"water": water, "r134a": r134a, "ethanol": ethanol, "tmp": tmp}

# Creates polynomial with Tc_o as the variable
def Tc_o_func(T, Th_i, Tc_i, Th_o, mh, mc, hf, cf):
    '''
    :param Th_i: input Temp for hot fluid
    :param Tc_i: input Temp for cold fluid
    :param Th_o: output Temp for hot fluid
    :param mh: mass flow rate of hot fluid
    :param mc: mass flow rate of cold fluid
    :param hf: str of hot fluid species
    :param cf: str of cold fluid species
    :param T: Temperature
    :return: a polynomial with Tc_o as the variable
    '''

    Th_avg = (Th_i + Th_o)/2.0  # average Temp of the hot fluid
    Cp_h = subst[hf]["A"] + subst[hf]["B"] * Th_avg + subst[hf]["C"] * Th_avg ** 2.0 \
           + subst[hf]["D"] * Th_avg ** 3.0 + subst[hf]["E"] * Th_avg ** 4.0  # heat capacity of hot fluid
    beta = subst[cf]["MW"] * mh * Cp_h * (Th_i - Th_o) / (mc*subst[hf]["MW"])
    a = subst[cf]["A"] + (1.0/4.0)*(-subst[cf]["C"]*Tc_i**2.0 - subst[cf]["D"] * Tc_i**3.0) \
        - 3.0*subst[cf]["E"]*Tc_i**4.0/16.0
    b = subst[cf]["B"]/2.0 + subst[cf]["C"]*Tc_i/4.0 - subst[cf]["E"]*Tc_i**3.0/8.0
    c = subst[cf]["C"]/4.0 + subst[cf]["D"]*Tc_i/4.0 + subst[cf]["E"]*Tc_i**2.0/8.0
    d = subst[cf]["D"]/8.0 + 3.0*subst[cf]["E"]*Tc_i/16.0
    e = subst[cf]["E"]/16.0
    alpha = -subst[cf]["A"]*Tc_i - subst[cf]["B"]*Tc_i**2.0/2.0 - subst[cf]["C"]*Tc_i**3.0/4.0\
            - subst[cf]["D"]*Tc_i**4.0/8.0 - subst[cf]["E"]*Tc_i**5.0/16.0
    y = a*T + b*T**2.0 + c*T**3.0 + d*T**4.0 + e*T**5.0 + alpha - beta
    return y

#creates polynomial with Th_o as the variable
def Th_o_func(T, Tc_i, Th_i, Tc_o, mh, mc, cf, hf):
    '''
    :param Tc_i: input Temp for cold fluid
    :param Th_i: input Temp for hot fluid
    :param Tc_o: output Temp for cold fluid
    :param mh: mass flow rate of hot fluid
    :param mc: mass flow rate of cold fluid
    :param cf: str of cold fluid species
    :param hf: str of hot fluid species
    :param T: temperature
    :return: a polynomial with Th_o as the variable
    '''

    Th_avg = (Tc_i + Tc_o) / 2.0  # average Temp of the cold fluid
    Cp_c = subst[cf]["A"] + subst[cf]["B"] * Th_avg + subst[cf]["C"] * Th_avg ** 2.0 \
           + subst[cf]["D"] * Th_avg ** 3.0 + subst[cf]["E"] * Th_avg ** 4.0  # heat capacity of cold fluid
    beta = subst[hf]["MW"] * mc * Cp_c * (Tc_o - Tc_i) / (mh * subst[cf]["MW"])
    a = -subst[hf]["A"] + (1.0 / 4.0) * (subst[hf]["C"] * Th_i ** 2.0 + subst[hf]["D"] * Th_i ** 3.0) \
        + 3.0 * subst[hf]["E"] * Th_i ** 4.0 / 16.0
    b = -subst[hf]["B"] / 2.0 - subst[hf]["C"] * Th_i / 4.0 + subst[hf]["E"] * Th_i ** 3.0 / 8.0
    c = -subst[hf]["C"] / 4.0 - subst[hf]["D"] * Th_i / 4.0 - subst[hf]["E"] * Th_i ** 2.0 / 8.0
    d = -subst[hf]["D"] / 8.0 - 3.0 * subst[hf]["E"] * Th_i / 16.0
    e = -subst[hf]["E"] / 16.0
    alpha = subst[hf]["A"] * Th_i + subst[hf]["B"] * Th_i ** 2.0 / 2.0 + subst[hf]["C"] * Th_i ** 3.0 / 4.0 +\
            subst[hf]["D"] * Th_i**4.0 / 8.0 + subst[hf]["E"] * Th_i ** 5.0 / 16.0
    y = a*T + b*T**2.0 + c*T**3.0 + d*T**4.0 + e*T**5.0 + alpha - beta
    return y


# Asks user to pick hot fluid species
hotf = input("Enter Hotter fluid species: ")

while True:
    Thi, Unit1 = input("Enter Hot Fluid inlet Temp & Units (F or K): ").split()
    if Unit1.upper() != "F" and Unit1.upper() != "K":
        print("Invalid Units.")
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
    if Thi < subst[hotf]["Mp"]:
        print("Invalid Temp. Solid substance.")
        continue
    elif Thi > subst[hotf]["Bp"]:
        print("Invalid Temp. Gaseous substance.")
        continue
    else:
        break

# Asks user to pick cold fluid species
coldf = input("Enter Colder fluid species: ")

while True:
    Tci, Unit2 = input("Enter Cold Fluid inlet Temp & Units (F or K): ").split()
    if Unit2.upper() != "F" and Unit2.upper() != "K":
        print("Invalid Units.")
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
    if Tci < subst[coldf]["Mp"]:
        print("Invalid Temp. Solid substance.")
        continue
    elif Tci > subst[coldf]["Bp"]:
        print("Invalid Temp. Gaseous substance.")
        continue
    else:
        break

# Asks user which outlet temp, Tc_o or Th_o, to specify
while True:
    To = input("Which outlet Temp to specify (Tc_o or Th_o)?: ")
    if To.lower() == "tc_o":
        while True:
            Tco, Unit3 = input("Enter Cold Fluid outlet temp & Units (F or K): ").split()
            if Unit3.upper() != "F" and Unit3.upper() != "K":
                print("Invalid Units.")
                continue
            try:
                Tco = float(Tco)
            except ValueError:
                print("Invalid Temp Entry")
                continue
            if Unit3.upper() == "F":
                Tco = tempSI(Tco)
            elif Unit3.upper() == "K":
                Tco = Tco
            if Tco < subst[coldf]["Mp"]:
                print("Invalid Temp. Solid substance.")
                continue
            elif Tco > subst[coldf]["Bp"]:
                print("Invalid Temp. Gaseous substance.")
                continue
            else:
                break
    elif To.lower() == "th_o":
        while True:
            Tho, Unit3 = input("Enter Hot Fluid outlet temp & Units (F or K): ").split()
            if Unit3.upper() != "F" and Unit3.upper() != "K":
                print("Invalid Units.")
                continue
            try:
                Tho = float(Tho)
            except ValueError:
                print("Invalid Temp Entry")
                continue
            if Unit3.upper() == "F":
                Tho = tempSI(Tho)
            elif Unit3.upper() == "K":
                Tho = Tho
            if Tho < subst[hotf]["Mp"]:
                print("Invalid Temp. Solid substance.")
                continue
            elif Tho > subst[hotf]["Bp"]:
                print("Invalid Temp. Gaseous substance.")
                continue
            else:
                break
    if To.lower() != "tc_o" and To.lower() != "th_o":
        print("Invalid Entry.")
        continue
    else:
        break

# Asks user for Mass flow rate of hot fluid
while True:
    mh, Unit4 = input("Enter Hot Fluid Mass Flow Rate & Units (lbm/s or kg/s): ").split()
    if Unit4.lower() != "lbm/s" and Unit4.lower() != "kg/s":
        print("Invalid Units.")
        continue
    try:
        mh = float(mh)
    except ValueError:
        print("Invalid Mass Flow Rate Entry.")
        continue
    if Unit4.lower() == "lbm/s":
        mh = mflowSI(mh)
    elif Unit4.lower() == "kg/s":
        mh = mh
    if mh < 0:
        raise SystemExit("BRUHHHH! Negative value for flow rate. PROGRAM TERMINATED.")
    else:
        break

# Asks user for Mass flow rate of Cold fluid
while True:
    mc, Unit4 = input("Enter Cold Fluid Mass Flow Rate & Units (lbm/s or kg/s): ").split()
    if Unit4.lower() != "lbm/s" and Unit4.lower() != "kg/s":
        print("Invalid Units.")
        continue
    try:
        mc = float(mc)
    except ValueError:
        print("Invalid Mass Flow Rate Entry.")
        continue
    if Unit4.lower() == "lbm/s":
        mc = mflowSI(mc)
    elif Unit4.lower() == "kg/s":
        mc = mc
    if mc < 0:
        raise SystemExit("BRUHHHH! Negative value for flow rate. PROGRAM TERMINATED.")
    else:
        break

# Asks user for Overall Heat Transfer Coefficient
while True:
    U, Unit5 = input("Enter Overall Heat Transfer Coefficient (U) & Units [Enter 1 for Btu/(hr-ft^2-F)."
                     " Enter 2 for J/(s-m^2-K)]: ").split()
    if Unit5 != "1" and Unit5!= "2":
        print("Invalid Units.")
        continue
    try:
        U = float(U)
    except ValueError:
        print("Invalid Overall Heat Transfer Coefficient.")
        continue
    if Unit5 == "1":
        U = U_SI(U)
    elif Unit5 == "2":
        U = U
    if U < 0:
        raise SystemExit("BRUHHHH! Negative value for flow rate. PROGRAM TERMINATED.")
    else:
        break

if To.lower() == "th_o":
    x1 = np.linspace(0, 1000, 2000)
    y1 = Tc_o_func(x1, Thi, Tci, Tho, mh, mc, hotf, coldf)
    y2 = np.zeros(len(x1))
    plt.plot(x1, y1, "b")
    plt.xlabel("Tc_o")
    plt.ylabel("y = f(Tc,o)")
    plt.plot(x1, y2, 'g')
    print("FROM PLOT, find guess value for root. Exit plot to continue with program.")
    plt.show()

elif To.lower() == "tc_o":
    x1 = np.linspace(0, 1000, 2000)
    y1 = Th_o_func(x1, Tci, Thi, Tco, mh, mc, coldf, hotf)
    y2 = np.zeros(len(x1))
    plt.plot(x1, y1, "b")
    plt.xlabel("Th,o")
    plt.ylabel("y = f(Th,o)")
    plt.plot(x1, y2, "g")
    print("FROM PLOT, find guess value for root. Exit plot to continue with program.")
    plt.show()

# Ask user for guess value for root
if To.lower() == "th_o":
    while True:
        guess = input("Enter guess value for Tc_o:")
        try:
            guess = float(guess)
        except ValueError:
            print("Invalid Guess")
            continue
        else:
            break

elif To.lower() == "tc_o":
    while True:
        guess = input("Enter guess value for Th_o:")
        try:
            guess = float(guess)
        except ValueError:
            print("Invalid Guess")
            continue
        else:
            break

# Solve for Tco if Th_o was specified or Tho if Tc_o was specified
if To.lower() == "th_o":
    Tco = fsolve(Tc_o_func, guess,args=(Thi, Tci, Tho, mh, mc, hotf, coldf))
elif To.lower() == "tc_o":
    Tho =fsolve(Th_o_func, guess, args=(Tci, Thi, Tco, mh, mc, coldf, hotf))

# Calculates Delta T lm
def get_delT_lm(Th_i, Th_o, Tc_i, Tc_o):
    delT1 = Th_i - Tc_o
    delT2 = Th_o - Tc_i
    return (delT2 - delT1)/np.log(delT2/delT1)

# Calculates F
def get_F(Th_i, Th_o, Tc_i, Tc_o):
    R = (Th_i - Th_o)/(Tc_o - Tc_i)
    P = (Tc_o - Tc_i)/(Th_i - Tc_i)
    a = np.sqrt(R**2 + 1)/(R - 1)
    b = np.log((1 - P) / (1 - P*R))
    c = np.log((2 - P*(R + 1 - np.sqrt(R**2 + 1))) / (2 - P*(R + 1 + np.sqrt(R**2 + 1))))
    return a*b/c

# Calculates A
def get_A(Th_i, Th_o, Tc_i, Tc_o, hf, u):
    Th_avg = (Th_i + Th_o) / 2.0  # average Temp of the hot fluid
    Cp_h = subst[hf]["A"] + subst[hf]["B"] * Th_avg + subst[hf]["C"] * Th_avg ** 2.0 + subst[hf]["D"] * Th_avg ** 3.0\
           + subst[hf]["E"] * Th_avg ** 4.0  # heat capacity of hot fluid
    F = get_F(Th_i, Th_o, Tc_i, Tc_o)
    DelT_lm = get_delT_lm(Th_i, Th_o, Tc_i, Tc_o)
    q = mh*Cp_h*(Th_i - Th_o)/(subst[hf]["MW"])
    A = q / (F * u * DelT_lm)
    return A

# Returns Deliverables
Area = get_A(Thi, Tho, Tci, Tco, hotf, U)

print("")
print("Surface area for heat exchanger: {SA} m^2 ".format(SA=Area))
print("Total cost for heat exchanger: $ {Cost} USD".format(Cost=1000*Area))

if To.lower() == "tc_o":
    print("Hot Fluid outlet Temp (Tho): {tho} K".format(tho=Tho))

elif To.lower() == "th_o":
    print("Cold Fluid outlet Temp (Tco): {tho} K".format(tho=Tco))
print("")

print("Your engineering skills are #OnFleek. Messi would be proud.")


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
