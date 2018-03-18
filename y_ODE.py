## Returns ODEs for y 
from utilities import *

def y_ODE(y,t,T,C,G,B,params) :


	rho = params["rho"]
	kv = params["kv"]
	#print t 
	

	
	dy1 = -3*rho*kv*G*(y[3]+y[7])
	dy2 = 0
	dy3 = G*y[1]
	dy4 =  2*G*y[2]
	dy5 = 3*G*y[3]
	dy6 = B
	dy7 = G*y[5]
	dy8 = 2*G*y[6]
	dy9 = 3*G*y[7]

	return [dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9]



