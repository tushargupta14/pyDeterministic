### Main File 
##--Input : Parameters 
##-- Executes all the functions 

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


def check_constraint(C,):



	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 


	if C <  Cs 
    ## Evaluat0e T from Cs expression
    p = [7.14*10^-6 -2.46*10^-3 (-6.29*10^-2+C)];
    disp('inside the Cs constraint');
    r = roots(p)
    if isreal(r(1)) ==0 
       C = Cs ;
       p = [7.14*10^-6 -2.46*10^-3 (-6.29*10^-2+C)];
       r = roots(p)
    end
    
    if min(r) > 0 
        
        T_new = min(r)+273
    else 
        T_new = max(r)+273
    end
    
    %Cm =  7.76 * 10^-2 + 2.46*10^-3 * (T_new-273) - 8.1*10^-6 * (T_new-273)^2 ;
    
    
    if C > Cm 
        
        p = [8.1*10^-6  -2.46*10^-3  (-7.76*10^-2+C)];
        disp('Inside the Cs Cm constraint');
        r = roots(p)
        
        t_min = min(r);
        if t_min > 0 
            T_new = t_min+273 
        else
            T_new = max(r) + 273
        end   
        %if 27 <= r(1) <= 52
            %T_new = r(1)+273
            
        %end
        %if  27<=r(2)<=52
            %T_new = r(2)+273
        %end
        %if 27<=r(1)<=52 && 27<=r(2)<=52 
            %T_new = min(r) +273
        
        %end
        
        if T_new == 0
            display('Error');
        end
    end
end
if T_new < 300 
    T_new = 300
end

if T_new > 325 
 T_new = 325

end

if T_new ==0 
 display('Here')
 T_new = T_computed;
end



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

	M = -10**-7
	tolerance = 10**-2

	num_iter = 2

	time_length = len(range(t0,tf+delta_t,delta_t))
	T_vec = np.ones(time_length)*323
	DH_vec = np.zeros((num_iter,time_length))



	iteration = 0


	print y0.shape


	while(iteration <= num_iter) :

		y_mat = np.zeros((time_length,9))
		z_mat = np.zeros((time_length,9))
		theta_mat = np.zeros((time_length,9))
		fi_mat = np.zeros((time_length,9))

		DelH_dy_mat = np.zeros((time_length,9))
		DelH_dz_mat = np.zeros((time_length,9))

		y_mat[0,:] = y0
		#print y_mat[0,0]
		z_mat[-1,:] = zf
		theta_mat[0,:] = theta0
		fi_mat[-1,:] = fi_f
	
		for t in range(t0,tf,delta_t) :
			#print t 
			t_horizon = np.linspace(t,t+delta_t,num = 10)
			#print t_horizon
			T = T_vec[t]
			C = y_mat[t,0]
			G = calG(T,C,parameters)
			B = calB(y_mat[t,:],T,parameters)

			y = odeint(y_ODE,y_mat[t,:],t_horizon,args = (T,C,G,B,parameters))
			
			y_mat[t+delta_t,:] = y[-1,:]

		for t in range(tf,t0,-delta_t):
			#print t 
			t_horizon = np.linspace(t,t-delta_t,num =10)
			#print t_horizon
			T = T_vec[t]
			C = y_mat[t,0]
			G = calG(T,C,parameters)
			z = odeint(z_ODE,z_mat[t,:],-t_horizon,args = (G,parameters,T,y_mat[t,:]))
			#print z[-1,0]			
			z_mat[t-delta_t,:] = z[-1,:]
		

		for t in range(t0,tf+delta_t,delta_t):

			##t_horizon = np.linspace(t,t+delta_t,num = 10)
			T = T_vec[t]
			G = calG(T,C,parameters)

			DelH_dy_mat[t,:] = DH_dy(y_mat[t,:],z_mat[t,:],G,T,parameters)
			DelH_dz_mat[t,:] = DH_dz(T,y_mat[t,:],parameters)
		
		## Theta forward integration
		for t in range(t0,tf,delta_t) :

			T = T_vec[t]
			t_horizon = np.linspace(t,t+delta_t,num = 10)
			theta = odeint(theta_ODE,theta_mat[t,:],t_horizon,args = (y_mat[t,:],T,parameters))

			theta_mat[t+delta_t,:] = theta[-1,:]


		for t in range(tf,t0,-delta_t):
			#print t 
			t_horizon = np.linspace(t,t-delta_t,num =10)
			#print t_horizon
			T = T_vec[t]
	
			fi = odeint(fi_ODE,fi_mat[t,:],-t_horizon,args = (y_mat[t,:],z_mat[t,:],theta_mat[t,:],T,parameters))
			#print z[-1,0]			
			fi_mat[t-delta_t,:] = fi[-1,:]
		

		#print fi_mat


		for t in range(t0,tf+delta_t,delta_t) :
			var_sum = 0 

			for i in range(9):
				var_sum = DelH_dy_mat[t,i]*theta_mat[t,i] + DelH_dz_mat[t,i]*fi_mat[t,i]

			DH_vec[iteration,t] = var_sum



		for t in range(t0,tf+delta_t,delta_t) :

			if DH_vec[iteration,t] < tolerance :

				T_vec[t] = check_constraint()




		iteration+=1
		break 

		## Fi backward Integration

	


	plt_z =  DH_vec[0,:]

	## Plotting function 

	t = np.linspace(t0,tf,num = 1801)
	plt.plot(t,plt_z)
	plt.show()	

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