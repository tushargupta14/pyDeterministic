### Utilities


import numpy as np 
import math

def Temp(t):

	##y = 45 - 20*(t/(240.0))**(0.33)
	##y = 45 - (20/240)*()
	x =  45.0 - (20.0)*(t/(240.0))
	return x

def SS(C,T):

	Cs = 5*10**(-5)*T**2 - 0.001*T + 0.0236
	if C > Cs :

		return (C-Cs)/Cs
	else :
		return 0

def calG(T,C,params):

	kg = params["kg"]

	g = params["g"]


	S = SS(C,T)
	#print S
	growth_rate = kg*(S**g)

	return growth_rate

def calB(y,T,params) :

	k_j1 = params["kj_1"]
	k_j2 = params["kj_2"]
	F = params["F"]
	Cc = params["Cc"]

	Cs = 5*(10**-5)*(T**2) - 0.001*T + 0.0236
	#S = y[0]/Cs
	#print -k_j2*((math.log(Cc/(Cs*(10**6)/F))**3)/(math.log(y[0]/Cs)**2))
	S = SS(y[0],T)

	if S!=0 :
		S = y[0]/Cs
	else :
		return 0 

	B = k_j1*S*np.exp(-k_j2*(math.log(Cc/(Cs*10**6/F))**3/math.log(S)**2))    ######
	#print B
	return B

def DG_dy(T,C,params):

	kg = params["kg"]

	g = params["g"]	

	Cs = 5*(10**-5)*(T**2) - 0.001*T + 0.0236
	
	S = SS(C,T)
	
	return kg*(S**(g-1))*g/Cs


def DB_dy(T,C,y,params):

	k_j1 = params["kj_1"]
	k_j2 = params["kj_2"]
	F = params["F"]
	Cc = params["Cc"]

	Cs = 5*(10**-5)*(T**2) - 0.001*T + 0.0236


	A = k_j1*(1/Cs)*np.exp(-k_j2*((math.log(Cc/(Cs*(10**6)/F))**3)/(math.log(y[0]/Cs)**2)))*60
	#print A
	B = k_j1*(y[0]/Cs)*np.exp(-k_j2*((math.log(Cc/(Cs*(10**6)/F))**3)/(math.log(y[0]/Cs)**2)))*60 
	#* (-k_j2)*((math.log(Cc/(Cs*(10**6)/F))**3)/(math.log(y[0]/Cs)**2))
	C = -2*(1/math.log(y[0]/Cs)**3)*(1/y[0])*B
	#print B
	return A+C

def DG_dt(theta,T,C,params):


	kg = params["kg"]

	g = params["g"]
	Cs = 5*(10**-5)*(T**2) - 0.001*T + 0.0236
	S  = SS(C,T)

	DCs_dT = 10*(10**-5)*(T) - 0.001 
	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (C - Cs))/Cs**2
	
	DelG_dT = kg*g*(S**(g-1))*exper

	return DelG_dT

def DB_dt(theta,y,T,params):


	Cs = 5*(10**-5)*(T**2) - 0.001*T + 0.0236
	DCs_dT =  10*(10**-5)*(T) - 0.001 	
	
	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (y[0] - Cs))/Cs**2

	k_j1 = params["kj_1"]
	k_j2 = params["kj_2"]
	F = params["F"]
	Cc = params["Cc"]
	C = y[0]
	S  = (C- Cs)/Cs

	return 0 

def DG_dydT(theta,T,C,params):


	Cs = Cs = 5*(10**-5)*(T**2) - 0.001*T + 0.0236
	kg = params["kg"]

	g = params["g"]

	DCs_dT = 10*(10**-5)*(T) - 0.001 	
	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (C - Cs))/Cs**2
	S  = (C- Cs)/Cs
	#A = kg*g*np.exp(-Eg/T)*S**(g-1)*(Eg/T**2)/Cs
	#B = kg*g*np.exp(-Eg/T)*(g-1)*S**(g-2)*exper/Cs
	#C = -1*kg*g*np.exp(-Eg/T)*S**(g-1)*DCs_dT/Cs**2

	return 0 

def DB_dydT(theta,y,T,params):


	Cs = Cs = 5*(10**-5)*(T**2) - 0.001*T + 0.0236
	DCs_dT = 10*(10**-5)*(T) - 0.001 	

	S  = (y[0]	- Cs)/Cs
	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (y[0] - Cs))/Cs**2


	return 0 
