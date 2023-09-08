function [p,t,rho]=atmosphere(z)
g=      9.80665;    %m/s^2
R=      287.04;     %J/kg/K
R_u=    8314.32;    %J/kmol/K
r_e=    6356766;    %m
M_0=    28.9644;
Z=      [0.001 11e3 2e4 32e3 47e3 51e3 71e3 84852];         %m
P=      [101325 22632 5474.8 868.01 110.9 66.94 3.69];      %Pa
T=      [288.15 216.65 216.65 228.65 270.65 270.65 214.65]; %K
S=      [-0.0065 0 0.0010 0.0028 0 -0.0028 -0.0020];        %K/m
for j=1:length(z)
    H=  r_e*z(j)/(r_e+z(j));
    a=  H./Z;
    i=  find(a>=1,1,'last');
    if isempty(i)
        i=   7; 
    end
    t(j)=   T(i)+S(i)*(H-Z(i));
    if i==2 || i==5
        p(j)=   P(i)*exp(-g*(H-Z(i))/R/T(i));
    else
        p(j)=   P(i)*(T(i)/t(j))^(g*M_0/R_u/S(i));
    end
end
rho=    p/R./t;