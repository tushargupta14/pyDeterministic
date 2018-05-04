### Main function for the LAM model


import numpy as np
from scipy.integrate import odeint
from y_ODE import *
from z_ODE import *
from DH_dy import *
from DH_dz import *
from theta_ODE import *
from fi_ODE import *
from utilities import *
import matplotlib.pyplot as plt
from check_constraint import *
import math



def model(parameters,delta_t = 1):

	t0 = 0 
	tf = 240 
	delta_t = 1


	kg = parameters["kg"] 
	g = parameters["g"] 
	kj_1 = parameters["kj_1"] 
	kj_2 = parameters["kj_2"] 
	rho = parameters["rho"] 
	kv = parameters["kv"]

		
	y0 = np.array([0.0732,0,0,0,0,0])

	y0 = y0.reshape(1,-1)

	zf = np.array([0,0,0,0,0,1])
	zf = zf.reshape(1,-1)

	theta0 = np.array([0,0,0,0,0,0])
	theta0 = theta0.reshape(1,-1)

	fi_f = np.array([0,0,0,0,0,0])
	fi_f = fi_f.reshape(1,-1)

	M = -10**-20
	tolerance = 10**-2

	num_iter = 1

	time_length = len(range(t0,tf+delta_t,delta_t))
	T_vec = np.array([Temp(t) for t in np.arange(t0,tf,delta_t)])
	DH_vec = np.zeros((num_iter,time_length))


	print time_length



	iteration = 0

	while(iteration < num_iter) :


		print iteration

		y_mat = np.zeros((time_length,6))
		z_mat = np.zeros((time_length,6))
		theta_mat = np.zeros((time_length,6))
		fi_mat = np.zeros((time_length,6))

		y_mat[0,:] = y0
		z_mat[0,:] = zf
		theta_mat[0,:] = theta0
		fi_mat[0,:] = fi_f

		DelH_dy_mat = np.zeros((time_length,6))
		DelH_dz_mat = np.zeros((time_length,6))

		print "y forward integration"
		
		for t in range(t0,tf,delta_t) :
			
			t_horizon = np.linspace(t,t+delta_t,num = 10)
			##print y_mat[t,:]
			T = T_vec[t]
			C = y_mat[t,0]
			G = calG(T,C,parameters)
			B = calB(y_mat[t,:],T,parameters)
			print B
			#print G
			dy_vec = y_ODE(y_mat[t,:],t,T,C,G,B,parameters)
			dy_vec = np.array([i*delta_t for i in dy_vec])
			#y = odeint(y_ODE,y_mat[t,:],t_horizon,args = (T,C,G,B,parameters))				
			y_mat[t+delta_t,:] = y_mat[t,:] + dy_vec
			#y = odeint(y_ODE,y_mat[t,:],t_horizon,args = (T,C,G,B,parameters))
			
			#y_mat[t+delta_t,:] = y[-1,:]
		"""
	 	print "Z backward ..."
		for t in range(t0,tf,delta_t):

			T = T_vec[t]
			C = y_mat[t,0]
			G = calG(T,C,parameters)

			t_horizon = np.linspace(t,t+delta_t,num = 10)
			z = odeint(z_ODE,z_mat[t,:],t_horizon,args = (G,parameters,T,y_mat[t,:]))
			z_mat[t+delta_t,:] = z[-1,:]
	

		print "DH .."
		for t in range(t0,tf,delta_t):

		
			T = T_vec[t]
			G = calG(T,C,parameters)

			DelH_dy_mat[t,:] = DH_dy(y_mat[t,:],z_mat[t,:],G,T,parameters)
			DelH_dz_mat[t,:] = DH_dz(T,y_mat[t,:],parameters)

		print " Theta forward integration.."
		for t in range(t0,tf,delta_t) :

			T = T_vec[t]
			t_horizon = np.linspace(t,t+delta_t,num = 10)
			theta = odeint(theta_ODE,theta_mat[t,:],t_horizon,args = (y_mat[t,:],T,parameters))

			theta_mat[t+delta_t,:] = theta[-1,:]

		
		print "Fi backward .."
		for t in range(t0,tf,delta_t):
			
			T = T_vec[t]
			t_horizon = np.linspace(t,t+delta_t,num = 10)


			fi = odeint(fi_ODE,fi_mat[t,:],t_horizon,args = (y_mat[t,:],z_mat[t,:],theta_mat[t,:],T,parameters))
			#print z[-1,0]			
			fi_mat[t+delta_t,:] = fi[-1,:]
	

		print "Derivative sums...."
		for t in range(t0,tf+delta_t,delta_t) :
			var_sum = 0

			for i in range(6):
				var_sum +=  DelH_dy_mat[t,i]*theta_mat[t,i] + DelH_dz_mat[t,i]*fi_mat[t,i]
				## + +
			DH_vec[iteration,t] = var_sum


	
		print "check constraints...."
		for t in range(t0,tf+delta_t,delta_t) :

			if abs(DH_vec[iteration,t]) > tolerance :

				T_vec[t] = check_constraint(T_vec[t],y_mat[t,0],DH_vec[iteration,t],M)
	
		"""

		iteration+=1
		t = np.linspace(t0,tf,num= 240)
		plt.plot(t,T_vec)
		plt.show()
		#print fi_mat


if __name__ == "__main__" :


	parameters = {}

	F = 150.14
	parameters["kg"] = np.exp(3.41)
	parameters["g"] = 1.48
	parameters["kj_1"] = np.exp(24.74)
	parameters["kj_2"] = 2.7*10**-2
	parameters["rho"] = 1568*10^3
	parameters["kv"] = 	math.pi/6
	parameters["Cc"] = 1568*10**3/F
	parameters["F"] = F
	model(parameters)