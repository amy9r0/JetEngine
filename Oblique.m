function [beta]=Oblique(M,theta,T)
format long
gam=    Gamma(T);
mu=     asin(1/M);
a=      tan(mu)^2;
b=      (gam*(1+a)+a-1)*tan(theta)/2;
c=      b+(a+1)*tan(theta);
d=      1-3*b.*c;
e=      27*a*b.^2+9*b.*c-2;
f=      sqrt(4*d.^3./e.^2-1);
beta=   atan2(((9*a*b+c)./(2*d)-e.*f./(6*b.*d).*tan(atan2(1,f)/3)),1);