## Returns ODE for the fi variable

from utilities import *
def fi_ODE(fi,t,y,T,params):

	Cs = 6.29 * (10**-2) + 2.46*(10**-3) * (T-273) - 7.14 * (10**-6) * (T-273)**2 
  
    S = (y[0] -Cs)/ Cs

    DelG_dely_delT = DG_dydT(theta,T,y[0],params)
    DelB_dely_delT = DelB_dydT(T,S,theta,y)
    DelG_dT = DG_dT(theta,T,y[0],params)
    
    A = fi(1)*(3*rho*kv*(y(end,4)+y(end,8))*DelG(S,T));
    B1 = z(end,1)*(3*rho*kv*(theta(end,4)+theta(end,8)))*DelG(S,T);
    C = z(end,1)*3*rho*kv*(y(end,4)+y(end,8))*DelG_dely_delT;
    
    D = -DelG(S,T)*(fi(3)*y(2)+z(3)*theta(2)+2*fi(4)*y(3)+2*z(4)*theta(3)+3*fi(5)*y(4)+3*z(5)*theta(4)+fi(7)*y(6)+z(7)*theta(6)+2*fi(8)*y(7)+2*z(8)*theta(7)+3*fi(9)*y(8)+3*z(9)*theta(8));
    E = -DelG_dely_delT*(z(end,3)*y(end,2)+2*z(end,4)*y(end,3)+3*z(end,5)*y(end,4)+z(end,7)*y(end,6)+2*z(end,8)*y(end,7)+3*z(end,9)*y(end,8));
    F = -fi(6)*DelB(S,y(end,4),y(end,8),T);
    G1 = -z(end,6)*DelB_dely_delT;
   
    Cs = 6.29 * 10^-2 + 2.46*10^-3 * T - 7.14 * 10^-6 * T^2 ;
    DCs_dT = 2.46*10^-3-14.28*10^-6*(T-273);
    exper = (Cs*(theta(1) - DCs_dT) - DCs_dT * (y(1) - Cs))/Cs^2;

    
    
    
    dtfi1 = A+B1+C+D+E+F+G1;
    dtfi2  = -(z(end,3)*DelG_dT+fi(3)*G);
    dtfi3 = -2*(fi(4)*G + z(end,4)*DelG_dT);
    dtfi4 = 3*rho*kv*DelG_dT - 3*z(end,5)*DelG_dT - 3*fi(5)*G;
    dtfi5 = -kb*(fi(6)*exp(-E_b/T)*S^b + z(end,6)*exp(-E_b/T)*S^b*E_b/T^2 + z(end,6)*exp(-E_b/T)*b*S^(b-1)*exper);
    dtfi6 = -(z(end,7)*DelG_dT+fi(7)*G);
    dtfi7 = -2*(fi(8)*G+z(end,8)*DelG_dT);
    dtfi8 = 3*rho*kv*z(1)*DelG_dT + 3*rho*kv*fi(1)*G -3*z(end,9)*DelG_dT-3*fi(9)*G;
    dtfi9 = -kb*(fi(6)*exp(-E_b/T)*S^b+z(end,6)*exp(-E_b/T)*S^b*E_b/T^2 + z(end,6)*exp(-E_b/T)*b*S^(b-1)*exper);
    
    
    
    
    dfi_dt = [dtfi1;dtfi2;dtfi3;dtfi4;dtfi5;dtfi6;dtfi7;dtfi8;dtfi9];