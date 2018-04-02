### Return DH_dz
### Return DH_dz
from utilities import * 


def DH_dz(T,y,params) :
	
	C = y[0]
	G = calG(T,C,params)
	B = calB(y,T,params)

	kv = params["kv"]
   	rho = params["rho"]
	dh_dz1 = -3*rho*kv*G*(y[3]+y[7])
	dh_dz2 = 0
	dh_dz3 = G*y[1]
	dh_dz4 = 2*G*y[2]
	dh_dz5 = 3*G*y[3] 
	dh_dz6 = B
	dh_dz7 = y[5]*G
	dh_dz8 = 2*G*y[6]
	dh_dz9 = 3*G*y[7]

	dh_dz =  [dh_dz1,dh_dz2,dh_dz3,dh_dz4,dh_dz5,dh_dz6,dh_dz7,dh_dz8,dh_dz9]

	return [i*10**-6 for i in dh_dz]


