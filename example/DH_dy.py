### Calculates DH_dy
from utilities import * 
def DH_dy(y,z,G,T,params) :

	C = y[0]
	B = calB(y,T,params)
	DelG = DG_dy(T,C,params)
   	DelB = DB_dy(T,C,y,params)

   	kv = params["kv"]
   	rho = params["rho"]

	dh_dy1 = DelG*(-z[0]*3*rho*kv*(y[3])+z[2]*y[1]+2*z[3]*y[2]+3*z[4]*y[3]+4*z[5]*y[4]) + (z[1]*DelB)
	dh_dy2 = z[2]*G
	dh_dy3 = z[3]*2*G
	dh_dy4 = z[4]*3*G - 3*rho*kv*G*z[0]
	dh_dy5 = 4*z[5]*G
	dh_dy6 = 0 


	DH_dy =  [dh_dy1,dh_dy2,dh_dy3,dh_dy4,dh_dy5,dh_dy6]

	return [i for i in DH_dy]
