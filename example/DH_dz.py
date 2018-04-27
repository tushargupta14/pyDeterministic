from utilities import * 


def DH_dz(T,y,params) :
	
	C = y[0]
	G = calG(T,C,params)
	B = calB(y,T,params)

	kv = params["kv"]
   	rho = params["rho"]
	dh_dz1 = -3*rho*kv*G*(y[3])
	dh_dz2 = B
	dh_dz3 = G*y[1]
	dh_dz4 = 2*G*y[2]
	dh_dz5 = 3*G*y[3] 
	dh_dz6 = 4*G*y[4]

	dh_dz =  [dh_dz1,dh_dz2,dh_dz3,dh_dz4,dh_dz5,dh_dz6]

	return [i*10**-6 for i in dh_dz]


