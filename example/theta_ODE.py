from utilities import *


def theta_ODE(theta,t,y,T,params):
   kv = params["kv"]
   rho = params["rho"]
   
   C = y[0]
   G = calG(T,C,params)
   DelG_dT = DG_dt(theta,T,C,params)
   DelB_dT = DB_dt(theta,y,T,params)


   dthe1 = -3*rho*kv*(theta[3])*G -3*rho*kv*(y[3])*DelG_dT
   dthe2 = DelB_dT
   dthe3 = theta[1]*G+y[1]*DelG_dT
   dthe4 = 2*y[2]*DelG_dT + 2*theta[2]*G
   dthe5 = 3*DelG_dT*y[3]+3*theta[3]*G
   dthe6 = 4*DelG_dT*y[3]+4*theta[3]*G

   dtheta =  [dthe1,dthe2,dthe3,dthe4,dthe5,dthe6]

   return [i for i in dtheta]
   #return dtheta