def y_ODE(y,t,T,C,G,B,params) :


	rho = params["rho"]
	kv = params["kv"]
	#print t



	dy0 = -3*rho*kv*G*(y[3])/10**6
	dy1 = B
	dy2 = G*y[1]
	dy3 = 2*G*y[2]
	dy4 = 3*G*y[3]
	dy5 = 4*G*y[4]
	return [dy0,dy1,dy2,dy3,dy4,dy5]
