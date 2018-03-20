 ## Returns ODE for theta
from utilities import *


def theta_ODE(theta,t,y,T,params):
   kv = params["kv"]
   rho = params["rho"]
   
   C = y[0]
   G = calG(T,C,params)
   DelG_dT = DG_dt(theta,T,C,params)
   DelB_dT = DB_dt(theta,y,T,params)


   dthe1 = -3*rho*kv*(theta[3]+theta[7])*G -3*rho*kv*(y[3]+y[7])*DelG_dT;
   dthe2 = 0
   dthe3 = theta[1]*G+y[1]*DelG_dT;
   dthe4 = 2*y[2]*DelG_dT + 2*theta[2]*G;
   dthe5 = 3*DelG_dT*y[3]+3*theta[3]*G;
   dthe6 = DelB_dT;
   dthe7 = theta[5]*G+DelG_dT*y[5];
   dthe8 = 2*G*theta[6]+ 2*DelG_dT*y[6];
   dthe9 = 3*G*theta[7] + 3*DelG_dT*y[7];
   
   return [dthe1,dthe2,dthe3,dthe4,dthe5,dthe6,dthe7,dthe8,dthe9]

