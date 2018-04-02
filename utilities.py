## Utility Functions
import numpy as np 
import math
def calG(T,C,params):

	kg = params["kg"]
	Eg = params["Eg"]

	g = params["g"]
	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	S  = (C- Cs)/Cs
	
	#print S**g, S, g 
	growth_rate = kg*np.exp(-Eg/T)*S**g
	#print growth_rate
	return growth_rate


def calB(y,T,params) :

	kb = params["kb"]
	Eb = params["Eb"]
	b = params["b"]
	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	S  = (y[0]- Cs)/Cs
	B = (kb*np.exp(-Eb/T))*(S**b)*(y[4]+y[8])    ######
	
	return B


def DG_dy(T,C,params):

	kg = params["kg"]
	Eg = params["Eg"]

	g = params["g"]	

	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	S  = (C- Cs)/Cs
	
	return kg*g*np.exp(-Eg/T)*(S**(g-1))/Cs


def DB_dy(T,C,y,params):

	kb = params["kb"]
	Eb = params["Eb"]

	b = params["b"]

	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	S  = (C- Cs)/Cs
	return kb*b*np.exp(-Eb/T)*(S**(b-1))*(y[3]+y[7])/Cs    ###

def DG_dt(theta,T,C,params):


	kg = params["kg"]
	Eg = params["Eg"]

	g = params["g"]
	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	S  = (C- Cs)/Cs

	DCs_dT = 2.46*10**-3 -14.28*(10**-6)*(T-273)
	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (C - Cs))/Cs**2
	A = kg*S**g*np.exp(-Eg/T)*Eg/(T**2)
	DelG_dT = A + kg*np.exp(-Eg/T)*g*(S**(g-1))*exper

	return DelG_dT

def DB_dt(theta,y,T,params):


	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	DCs_dT = 2.46*10**-3 -14.28*10**-6*(T-273)	
	
	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (y[0] - Cs))/Cs**2

	kb = params["kb"]
	Eb = params["Eb"]
	C = y[0]
	S  = (C- Cs)/Cs


	b = params["b"]
	A = kb*np.exp(-Eb/T)*Eb*S**b*(y[4]+y[8])/T**2 #######
	B = A + kb*np.exp(-Eb/T)*b*S**(b-1)*exper*(y[4]+y[8]) ####
	DelB_dT = kb*np.exp(-Eb/T)*S**b*(theta[4]+theta[8])+B ####

	return DelB_dT

def DG_dT(theta,T,C,params):



	kg = params["kg"]
	Eg = params["Eg"]

	g = params["g"]

	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 

	DCs_dT = 2.46*10**-3 -14.28*10**-6*(T-273)	

	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (C - Cs))/Cs**2
	S  = (C- Cs)/Cs
	A = kg*S**g*np.exp(-Eg/T)*Eg/(T**2)

	DelG_dT = A + kg*np.exp(-Eg/T)*g*(S**(g-1))*exper

	return DelG_dT


def DG_dydT(theta,T,C,params):


	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	kg = params["kg"]
	Eg = params["Eg"]

	g = params["g"]

	DCs_dT = 2.46*10**-3 -14.28*10**-6*(T-273)
	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (C - Cs))/Cs**2
	S  = (C- Cs)/Cs
	A = kg*g*np.exp(-Eg/T)*S**(g-1)*(Eg/T**2)/Cs
	B = kg*g*np.exp(-Eg/T)*(g-1)*S**(g-2)*exper/Cs
	C = -1*kg*g*np.exp(-Eg/T)*S**(g-1)*DCs_dT/Cs**2

	return  A+B+C

def DB_dydT(theta,y,T,params):



	kb = params["kb"]
	Eb = params["Eb"]

	b = params["b"]
	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 

	DCs_dT = 2.46*10**-3 -14.28*10**-6*(T-273)

	S  = (y[0]	- Cs)/Cs
	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (y[0] - Cs))/Cs**2

	A = kb*b*np.exp(-Eb/T)*(Eb/T**2)*S**(b-1)*(y[3]+y[7])/Cs
	B = kb*b*np.exp(-Eb/T)*(b-1)*S**(b-2)*(y[3]+y[7])*exper
	C = -kb*b*np.exp(-Eb/T)*S**(b-1)*DCs_dT*(y[3]+y[7])/Cs**2            
	D = kb*b*np.exp(-Eb/T)*S**(b-1)*(theta[3]+theta[7])               

	return A+B+C+D


def solve_quadratic(a,b,c):


	d = math.sqrt(b**2 - 4*a*c)

	
	r1 = (-b + d)/(2*a)

	r2= (-b -d )/(2*a) 


	return r1,r2


