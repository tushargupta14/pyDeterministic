## The script contains the deterministic model 

## --  Takes input the delta_t , model parameters

## --  Output T_vec and y_mat




### Consists of odeint for all forward integeration
## This is the third testing script



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


   
def model(parameters,delta_t = 1,):


	t0 = 0
	tf = 1800  ## batch_time

	kg = parameters["kg"]
	Eg = parameters["Eg"]
	g = parameters["g"]
	kb = parameters["kb"]
	Eb = parameters["Eb"]
	b = parameters["b"]
	#rho = parameters["rho"]


	y0 = np.array([0.1743,66.66,1.83*10**4,5.05*10**6,1.93*10**9,0.867,0,0,0])
	y0 = y0.reshape(1,-1)
	zf = np.array([0,0,0,0,1,0,0,0,-1])
	zf = zf.reshape(1,-1)
	theta0 = np.array([0,0,0,0,0,0,0,0,0])
	theta0 = theta0.reshape(1,-1)
	fi_f = np.array([0,0,0,0,0,0,0,0,0])
	fi_f = fi_f.reshape(1,-1)

	M = -10**-6
	tolerance = 10**-2

	num_iter = 2

	time_length = len(range(t0,tf+delta_t,delta_t))
	T_vec = np.ones(time_length)*323
	DH_vec = np.zeros((num_iter,time_length))



	iteration = 0


	print y0.shape


	while(iteration < num_iter) :


		print iteration
		y_mat = np.zeros((time_length,9))
		z_mat = np.zeros((time_length,9))
		theta_mat = np.zeros((time_length,9))
		fi_mat = np.zeros((time_length,9))

		DelH_dy_mat = np.zeros((time_length,9))
		DelH_dz_mat = np.zeros((time_length,9))

		y_mat[0,:] = y0
		#print y_mat[0,0]
		z_mat[0,:] = zf
		theta_mat[0,:] = theta0
		fi_mat[0,:] = fi_f
		

		print "y forward integration"
		for t in range(t0,tf,delta_t) :
			
			t_horizon = np.linspace(t,t+delta_t,num = 10)
			T = T_vec[t]
			C = y_mat[t,0]
			G = calG(T,C,parameters)
			B = calB(y_mat[t,:],T,parameters)


			y = odeint(y_ODE,y_mat[t,:],t_horizon,args = (T,C,G,B,parameters))
			
			y_mat[t+delta_t,:] = y[-1,:]

		print "Z backward ..."
		for t in range(t0,tf,delta_t):

			T = T_vec[t]
			C = y_mat[t,0]
			G = calG(T,C,parameters)

			t_horizon = np.linspace(t,t+delta_t,num = 10)
			z = odeint(z_ODE,z_mat[t,:],t_horizon,args = (G,parameters,T,y_mat[t,:]))
			z_mat[t+delta_t,:] = z[-1,:]
	
		print "DH .."
		for t in range(t0,tf+delta_t,delta_t):

		
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

			for i in range(9):
				var_sum +=  DelH_dy_mat[t,i]*theta_mat[t,i] + DelH_dz_mat[t,i]*fi_mat[t,i] 
				## + +  
			DH_vec[iteration,t] = var_sum

		print "check constraints...."
		for t in range(t0,tf+delta_t,delta_t) :

			if abs(DH_vec[iteration,t]) > tolerance :

				#print "Here"
				T_vec[t] = check_constraint(T_vec[t],y_mat[t,0],DH_vec[iteration,t],M)


		plt_1 =  DH_vec[iteration,:]

	
		## Plotting function 
		"""
		t = np.linspace(t0,tf,num = 1801)
		plt.figure(0)
		plt.plot(t,plt_1,'b')
		#plt.show()
		plt.figure(1)
		plt.plot(t,T_vec)
		plt.figure(2)
		plt.plot(t,y_mat[:,4]-y_mat[:,8])
		plt.show()
		plt.cla()

		"""

		iteration+=1
				 


	return T_vec,y_mat
	#print T_vec
	"""
	plt_1 =  T_vec

	## Plotting function 

	t = np.linspace(t0,tf,num = 1801)
	#plt.plot(t,plt_1,'r')
	plt.figure(0)
	plt.plot(t,plt_1,'b')
	#plt.show()
	#plt.figure(1)
	#plt.plot(t,T_vec)
	plt.show()
	"""
if __name__ == "__main__" :

	parameters = {}
	parameters["kg"] = 1.44*10**8
	parameters["Eg"] = 4859 
	parameters["g"] = 1.5
	parameters["kb"] = 285
	parameters["Eb"] = 7517
	parameters["b"] = 1.45
	parameters["rho"] = 2.66*10**-12;
	parameters["kv"] = 0.54
	model(parameters)
