## Utility Functions
import numpy as np 

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
	B = (kb*np.exp(-Eb/T))*(S**b)*(y[4]+y[8])
	
	return B


def DG_dy(T,C,params):

	kg = params["kg"]
	Eg = params["Eg"]

	g = params["g"]	

	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	S  = (C- Cs)/Cs
	
	return kg*g*np.exp(-Eg/T)*(S**(g-1))/Cs


def DB_dy(T,C,y4,y8,params):

	kb = params["kb"]
	Eb = params["Eb"]

	b = params["b"]

	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	S  = (C- Cs)/Cs
	return kb*b*np.exp(-Eb/T)*(S**(b-1))*(y4+y8)/Cs