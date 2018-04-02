## Returns ODE for the fi variable

from utilities import *
import numpy as np 
def fi_ODE(fi,t,y,z,theta,T,params):

	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
	b = params["b"]
	kb = params["kb"]
	Eb = params["Eb"]
	kv = params["kv"]
	g = params["g"]
	kg = params["kg"]
	Eg = params["Eg"]
   	rho = params["rho"]

	S = (y[0]-Cs)/Cs
	G = calG(T,y[0],params)

	DelG_dely_delT = DG_dydT(theta,T,y[0],params)
	DelB_dely_delT = DB_dydT(theta,y,T,params)
	DelG_dT = DG_dT(theta,T,y[0],params)
	C = y[0]
	DelG = DG_dy(T,C,params)
	DelB = DB_dy(T,C,y,params)

	A = fi[0]*(3*rho*kv*(y[3]+y[7])*DelG)
	B1 = z[0]*(3*rho*kv*(theta[3]+theta[7]))*DelG
	C = z[0]*3*rho*kv*(y[3]+y[7])*DelG_dely_delT

	D = -DelG*(fi[2]*y[1]+z[2]*theta[1]+2*fi[3]*y[2]+2*z[3]*theta[2]+3*fi[4]*y[3]+3*z[4]*theta[3]+fi[6]*y[5]+z[6]*theta[5]+2*fi[7]*y[6]+2*z[7]*theta[6]+3*fi[8]*y[7]+3*z[8]*theta[7])
	E = -DelG_dely_delT*(z[2]*y[1]+2*z[3]*y[2]+3*z[4]*y[3]+z[6]*y[5]+2*z[7]*y[6]+3*z[8]*y[7])
	F = -fi[5]*DelB
	G1 = -z[5]*DelB_dely_delT

	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 

	DCs_dT = 2.46*10**-3 -14.28*10**-6*(T-273)
	exper = (Cs*(theta[0] - DCs_dT) - DCs_dT * (y[0] - Cs))/(Cs**2)


    
		   
	dtfi1 = A+B1+C+D+E+F+G1
	dtfi2  = -(z[2]*DelG_dT+fi[2]*G)
	dtfi3 = -2*(fi[3]*G + z[3]*DelG_dT)
	dtfi4 = 3*rho*kv*DelG_dT - 3*z[4]*DelG_dT - 3*fi[4]*G
	dtfi5 = -kb*(fi[5]*np.exp(-Eb/T)*S**b + z[5]*np.exp(-Eb/T)*(S**b)*Eb/(T**2) + z[5]*np.exp(-Eb/T)*b*S**(b-1)*exper)
	dtfi6 = -(z[6]*DelG_dT+fi[6]*G)
	dtfi7 = -2*(fi[7]*G+z[7]*DelG_dT)
	dtfi8 = 3*rho*kv*z[0]*DelG_dT + 3*rho*kv*fi[0]*G -3*z[8]*DelG_dT-3*fi[8]*G
	dtfi9 = -kb*(fi[5]*np.exp(-Eb/T)*(S**b)+z[5]*np.exp(-Eb/T)*S**b*Eb/(T**2) + z[5]*np.exp(-Eb/T)*b*S**(b-1)*exper)




	return [dtfi1,dtfi2,dtfi3,dtfi4,dtfi5,dtfi6,dtfi7,dtfi8,dtfi9]