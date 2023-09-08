function [L_1,L_2,L_3,vol,y_i]=Coordinates(th_1,th_2,th_3,b_1,b_2,b_3,b_4,L_i,L_b,plots,th_n,M)
%th-    wedge angle
%b-     mach angle
%p-     all points at change in wedge angle and point of cowl
p_c=    [0,-tan(th_1+b_1)*L_i];
m1=     -tan(th_1);
c1=     m1*L_i;
n1=     -tan(th_1+th_2+b_2);
d1=     p_c(2);
p_2=    [(-c1+d1)/(m1-n1),(d1*m1-c1*n1)/(m1-n1)];
f1=     @(x)(m1*x+c1);
m2=     -tan(th_1+th_2);
c2=     p_2(2)-m2*p_2(1);
n2=     -tan(th_1+th_2+th_3+b_3);
d2=     L_i*-tan(th_1+b_1);
p_3=    [(-c2+d2)/(m2-n2),(d2*m2-c2*n2)/(m2-n2)];
f2=     @(x)(m2*x+c2);
m3=     -tan(th_1+th_2+th_3);
c3=     p_3(2)-m3*p_3(1);
n3=     tan(b_4);
d3=     p_c(2);
p_4=    [(-c3+d3)/(m3-n3),(d3*m3-c3*n3)/(m3-n3)];
y_i=    p_4(2);
f3=     @(x)(m3*x+c3);
L_1=    norm([-L_i,0]-p_2);
L_2=    norm(p_2-p_3);
L_3=    norm(p_3-p_4);
vol=    0.5*p_2(1)*p_2(2)+0.5*(-p_2(2)-p_3(2))*(p_3(1)-p_2(1))+...
        0.5*(-p_3(2)-p_4(2))*(p_4(1)-p_3(1))-p_4(2)*L_b;
m4=     tan(th_n);
c4=     p_4(2)-m4*(p_4(1)+L_b);
f4=     @(x)(m4*x+c4);
if plots
    figure
    x=  linspace(-L_i,p_2(1),5);
    plot(x+L_i,f1(x),'k'), hold on
    x=  linspace(p_2(1),p_3(1),5);
    plot(x+L_i,f2(x),'k'), hold on
    x=  linspace(p_3(1),p_4(1),5);
    plot(x+L_i,f3(x),'k'), hold on
    plot(linspace(p_c(1),p_4(1),5)+L_i,p_c(2)*ones(1,5),'k')
    hold on
    plot(linspace(-L_i,-c4/m4,5)+L_i,zeros(1,5),'k')
    hold on, plot(p_2(1)+L_i,p_2(2),'.k')
    hold on, plot(p_3(1)+L_i,p_3(2),'.k')
    hold on, plot(linspace(p_4(1),p_4(1)+L_b,5)+L_i,zeros(1,5),'k')
    hold on, plot(linspace(p_4(1),p_4(1)+L_b,5)+L_i,ones(1,5)*p_4(2),'k')
    hold on, plot(linspace(p_4(1),p_4(1)+L_b,5)+L_i,ones(1,5)*p_c(2),'k')
    x=  linspace(p_4(1)+L_b,-c4/m4,5);
    hold on, plot(x+L_i,f4(x),'k')
    grid, axis equal
    xlabel('m'), ylabel('m')
    title({['Aircraft at M = ',num2str(M),'.']})
end