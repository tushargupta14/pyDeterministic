### Temperature constraint checking

from utilities import *

def check_constraint(T,C,DH_dt,M):






	T_computed = T + M*DH_dt

	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T_computed-273) - 7.14 * (10**-6) * (T_computed-273)**2 
	Cm =  7.76 * 10**-2 + 2.46*10**-3 * (T_computed-273) - 8.1*10**-6 * (T_computed-273)**2 
    

	T_new = 0


	if Cs < C :

		if C < Cm :

			T_new = T_computed

		

		else :


			a = 8.1*10**-6
			b = -2.46*10**-3
			c = (-7.76*10**-2+C)

			r1,r2 = solve_quadratic(a,b,c)
				
			

			if r2 > 0 :

				T_new = r2 + 273 
			else :
				T_new = r1+273	





	if C <  Cs :
    ## Evaluat0e T from Cs expression

	    a = 7.14*10**-6
	    b = -2.46*10**-3
	    c = (-6.29*10**-2+C)


	    r1,r2 = solve_quadratic(a,b,c)

	    if isinstance(r1,complex):

	    	print " 1"
	    	C = Cs
	    	c = (-6.29*10**-2+C)

	    	r1,r2 = solve_quadratic(a,b,c)
	 	

	    if r2 > 0 :

	    	T_new = r2 + 273 
	    else :
	    	T_new = r1+273



    
		if (C > Cm) :

			a = 8.1*10**-6
			b = -2.46*10**-3
			c = (-7.76*10**-2+C)

			r1,r2 = solve_quadratic(a,b,c)
			if isinstance(r1,complex):

				print "2"
				C = Cm
				c = (-7.76*10**-2+C)

				r1,r2 = solve_quadratic(a,b,c)
				

			if r2 > 0 :

				T_new = r2 + 273 
			else :
				T_new = r1+273	


	if T_new < 303:

		T_new = 303

	if T_new < 0 :

		print "egfefefefe"


	if T_new > 325 :

		T_new = 325 



	return T_new
 
    