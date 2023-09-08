function [I]=CompositeSimpson(a,b,y,n);
%a- lower bound of integral
%b- upper bound of integral
%y- function values at evenly spaced intervals
%n- number of evenly spaced intervals (even)
h=      (b-a)/n;
sum1=   0;
sum2=   0;
for i=1:n/2-1
    sum1=   sum1+y(2*i);
end
for j=1:n/2
    sum2=   sum2+y(2*j-1);
end
I=  h/3*(y(1)+2*sum1+4*sum2+y(end));