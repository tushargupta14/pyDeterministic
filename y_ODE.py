## Returns ODEs for y 


def y_ODE(y,T_vec,t,params) :


	rho = params["rho"]
	kv = params["kv"]

	dy1 = -3*rho*kv*G*(y[3]+y[7])
