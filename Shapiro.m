function [T0_e,M_sq_e,T_e,P_e]=Shapiro(x0,h,IC,tau_b,L,theta)
%x0-    initial burner value
%h-     equi-step value for range
%IC-    initial conditions,T0,M^2,T,P, column vector
%tau_b- total temperature ratio of burner (outlet/inlet)
%L-     length of burner
%_e-    exit burner properties
dT0=        @(x)(IC(1)*(tau_b-1)*(theta*L/(L^2+2*(theta-1)*L*x+(theta-1)*x^2)));
f1=         @(x,y)(y(1)*(tau_b-1)*(theta*L/(L^2+2*(theta-1)*L*x+(theta-1)*x^2)));%                  T0
f2=         @(x,y)(y(2)*(1+Gamma(y(3))*y(2))*(1+(Gamma(y(3))-1)/2*y(2))/(1-y(2))/y(1)*dT0(x));%     M^2
f3=         @(x,y)(y(3)*(1-Gamma(y(3))*y(2))*(1+(Gamma(y(3))-1)/2*y(2))/(1-y(2))/y(1)*dT0(x));%     T
f4=         @(x,y)(y(4)*-Gamma(y(3))*y(2)*(1+(Gamma(y(3))-1)/2*y(2))/(1-y(2))/y(1)*dT0(x));%        P
f=          {f1,f2,f3,f4};
[x,y]=      Runge_Kutta_4(x0,IC,h,L,f);
T0_e=       y(1,end);
M_sq_e=     y(2,end);
T_e=        y(3,end);
P_e=        y(4,end);