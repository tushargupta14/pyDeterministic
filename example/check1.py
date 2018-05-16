def check_constraint(T,C,DH_dt,M,t,e):



    T_computed = T + M*DH_dt

    Cs = 5*(10**-5)*(T**2) - 0.001*T + 0.0236
    #Cm =  7.76 * 10**-2 + 2.46*10**-3 * (T_computed-273) - 8.1*10**-6 * (T_computed-273)**2


    T_new = T

    if t <=25:
        try :
            T_new = 45.0 - (10.0)*(t/(25.0))**(e)
        except Exception as e :
            T_new = 45.0 - (10.0)*(t/(25.0))**(0.6)
    if t>25 :

        T_new = 34.0 - (2.5)*(t/(125.0))

    if t>=150 :
        T_new = 31 - (6.0)*((t-150)/(100.0))**(2)

    return T_new
