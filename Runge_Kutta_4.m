function [x,y]=Runge_Kutta_4(x0,y0,h,xf,f)
%x0-    initial range value
%y0-    initial function value(s), column vector
%h-     resolution of range values
%xf-    final range value
%f-     cell of functions
%x-     vector of evaluated range values
%y-     matrix of function values, each row corresponds to a function
x=      x0:h:xf;
y=      zeros(length(f),length(x));
y(:,1)= y0;
for i=2:length(x)
    for j=1:length(f)
        k_1(j,1)=   h*f{j}(x(i-1),y(:,i-1));
    end
    for m=1:length(f)
        r=          y(:,i-1)+0.5*k_1;
        k_2(m,1)=   h*f{m}(x(i-1)+h/2,r);
    end
    for n=1:length(f)
        s=          y(:,i-1)+0.5*k_2;
        k_3(n,1)=   h*f{n}(x(i-1)+h/2,s);
    end
    for p=1:length(f)
        u=          y(:,i-1)+k_3;
        k_4(p,1)=   h*f{p}(x(i-1)+h,u);
    end
    y(:,i)=   y(:,i-1)+(k_1+2*k_2+2*k_3+k_4)/6;
end