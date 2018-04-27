### z_ODE


from utilities import *
import numpy as np 
def z_ODE(z,t,G,params,T,y):


	rho = params["rho"]
	kv = params["kv"]



   	C = y[0]
   	Cs = 5*(10**-5)*(T**2) - 0.001*T + 0.0236


	S  = SS(C,T)

   	DelG = DG_dy(T,C,params)
   	DelB = DB_dy(T,C,y,params)


   	dz1 = z[0]*(3*rho*kv*(y[3])*DelG)

   	dz1 = dz1 - (z[2]*y[1]+2*z[3]*y[2]+3*z[4]*y[3]+4*z[5]*y[4])*DelG-DelB*z[1]  
   	dz2 = -z[2]*G
   	dz3 = -2*z[3]*G
	dz4 =  3*z[0]*rho*kv*G -3*z[4]*G
	dz5 = -4*z[5]*G
	dz6 = 0
	

	return [dz1,dz2,dz3,dz4,dz5,dz6]