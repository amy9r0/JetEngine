function [ve_v0,P_e]=Turbojet(T_inf,P_inf,h_PR,f,tau_b,M_inf,tau_c)
%T_0inf-    freestream total temperature
%T_inf-     freestream temperature
%T_3-       burner inlet temperature
R=          287;%                           J/kg/K
gam=        Gamma(T_inf);
c_p=        gam*R/(gam-1);
T_0inf=     T_inf*(1+(gam-1)/2*M_inf^2);
theta_inf=  T_0inf/T_inf;
theta_n=    f*h_PR/c_p/T_inf+theta_inf;
T_0_4=      tau_b*tau_c*T_0inf;
theta_t=    T_0_4/T_inf;
ve_v0=      sqrt(theta_n/(theta_n-1)*(1-theta_t/(theta_inf*tau_c)/(theta_t-theta_inf*(tau_c-1))));
T_e=        theta_n*T_inf;
M_e=        ve_v0*M_inf*sqrt(gam*T_inf/Gamma(T_e)/T_e);
P_e=        P_inf;