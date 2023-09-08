function [dT0_dx]=TotalTempChange(x,tau_b,T01,L,theta)
%x-         location along burner
%tau_b-     total temperature change across burner
%T01-       initial burner total temperature
%L-         total length of burner
%theta-     characteristic constant of burner
%dT0_dx-    change in total temperature with respect to burner location
dT0_dx= T01*(tau_b-1)*(theta*L/(L^2+2*(theta-1)*L*x+(theta-1)*x^2));