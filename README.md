# Chen-263-Heat-Exchanger-Project
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
