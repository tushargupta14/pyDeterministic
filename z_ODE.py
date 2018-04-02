## The ddifferential equations to solve z 
from utilities import *
import numpy as np 
def z_ODE(z,t,G,params,T,y):


	b = params["b"]
	kb = params["kb"]
	Eb = params["Eb"]
	kv = params["kv"]
	g = params["g"]
	kg = params["kg"]
	Eg = params["Eg"]
   	rho = params["rho"]
   	C = y[0]
   	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	S  = (C- Cs)/Cs
   	DelG = DG_dy(T,C,params)
   	DelB = DB_dy(T,C,y,params)


   	dz1 = z[0]*(3*rho*kv*(y[3]+y[7])*DelG)

   	dz1 = dz1 - (z[2]*y[1]+2*z[3]*y[2]+3*z[4]*y[3]+z[6]*y[5]+2*z[7]*y[6]+3*z[8]*y[7])*DelG-DelB*z[5]  
   	dz2 = -z[2]*G
   	dz3 = -2*z[3]*G
	dz4 =  3*z[0]*rho*kv*G -3*z[4]*G
	dz5 = -z[5]*kb*np.exp(-Eb/T)*S**b
	dz6 = -z[6]*G
	dz7 = -2*G*z[7]
	dz8 =  3*z[0]*rho*kv*G-3*z[8]*G
	dz9 = -z[5]*kb*np.exp(-Eb/T)*S**b

	return [dz1,dz2,dz3,dz4,dz5,dz6,dz7,dz8,dz9]