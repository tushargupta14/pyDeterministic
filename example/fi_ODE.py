## Returns ODE for the fi variable

from utilities import *
import numpy as np 
def fi_ODE(fi,t,y,z,theta,T,params):

	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 

	kv = params["kv"]
	g = params["g"]
	kg = params["kg"]

   	rho = params["rho"]

	S = (y[0]-Cs)/Cs
	G = calG(T,y[0],params)

	DelG_dely_delT = DG_dydT(theta,T,y[0],params)
	DelB_dely_delT = DB_dydT(theta,y,T,params)
	DelG_dT = DG_dt(theta,T,y[0],params)
	C = y[0]
	DelG = DG_dy(T,C,params)
	DelB = DB_dy(T,C,y,params)

	A = fi[0]*(3*rho*kv*(y[3])*DelG)
	B1 = z[0]*(3*rho*kv*(theta[3]))*DelG
	C = z[0]*3*rho*kv*(y[3])*DelG_dely_delT

	D = -DelG*(fi[2]*y[1]+z[2]*theta[1]+2*fi[3]*y[2]+2*z[3]*theta[2]+3*fi[4]*y[3]+3*z[4]*theta[3])
	E = -DelG_dely_delT*(z[2]*y[1]+2*z[3]*y[2]+3*z[4]*y[3])
	F = -fi[5]*DelB
	G1 = -z[5]*DelB_dely_delT

	Cs = 5*(10**-5)*(T**2) - 0.001*T + 0.0236
	DCs_dT =  10*(10**-5)*(T) - 0.001 	

	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (y[0] - Cs))/(Cs**2)

		   
	dtfi1 = A+B1+C+D+E+F+G1
	dtfi2  = -(z[2]*DelG_dT+fi[2]*G)
	dtfi3 = -2*(fi[3]*G + z[3]*DelG_dT)
	dtfi4 = 3*rho*kv*DelG_dT - 3*z[4]*DelG_dT - 3*fi[4]*G
	dtfi5 = 0 
	dtfi6 = -(z[5]*DelG_dT+fi[5]*G)



	return [dtfi1,dtfi2,dtfi3,dtfi4,dtfi5,dtfi6]