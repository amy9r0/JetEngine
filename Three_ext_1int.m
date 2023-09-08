function [P_e,T_e,M_e,beta_1,beta_2,beta_3,beta_4,P_2,P_3,P_4]...
            =Three_ext_1int(theta_1,theta_2,theta_3,theta_4,M,P,T)
%theta_1-   first shock wedge angle
%theta_2-   second shock wedge angle
%theta_3-   third shock wedge angle
%M-         freestream mach
%P-         freestream pressure
%T-         freestream temperature
%P_e-       exit pressure of inlet
%T_e-       exit temperature of inlet
%M_e-       exit mach of inlet

format long
gam=    Gamma(T);
%first shock (external) Mn_1=Mn_2
beta_1= Oblique(M,theta_1,T);
Mn_1=   M*sin(beta_1);
PR_1=   1+2*gam/(gam+1)*(Mn_1^2-1);
dR_1=   (gam+1)*Mn_1^2/((gam-1)*Mn_1^2+2);
TR_1=   PR_1/dR_1;
T_2=    TR_1*T;
gam=    Gamma(T_2);
Mn_2=   sqrt(((gam-1)*Mn_1^2+2)/(2*gam*Mn_1^2-gam+1));
M_2=    Mn_2/sin(beta_1-theta_1);
P_2=    PR_1*P;

%second shock (external) PR_1=PR_2 Mn_1=Mn_2=Mn_3 dR_1=dR_2 TR_1=TR_2
%beta_2=     asin(M/M_2*sin(beta_1));
beta_2=     Oblique(M_2,theta_2,T_2);
Mn_2=       M_2*sin(beta_2);
dR_2=       (gam+1)*Mn_2^2/((gam-1)*Mn_2^2+2);
PR_2=       1+2*gam/(gam+1)*(Mn_2^2-1);
TR_2=       PR_2/dR_2;
%theta_2=    atan2((2*cot(beta_2)*(Mn_2^2-1))/(M_2^2*(gam+cos(2*beta_2))+2),1);
T_3=        TR_2*T_2;
gam=        Gamma(T_3);
Mn_3=       sqrt(((gam-1)*Mn_2^2+2)/(2*gam*Mn_2^2-gam+1));
M_3=        Mn_3/sin(beta_2-theta_2);
P_3=        PR_2*P_2;

%third shock (external) PR_1=PR_3 Mn_1=Mn_3 dR_1=dR_3 TR_1=TR_3
%beta_3=     asin(M_2/M_3*sin(beta_2));
beta_3=     Oblique(M_3,theta_3,T_3);
Mn_3=       M_3*sin(beta_3);
dR_3=       (gam+1)*Mn_3^2/((gam-1)*Mn_3^2+2);
PR_3=       1+2*gam/(gam+1)*(Mn_3^2-1);
TR_3=       PR_3/dR_3;
%theta_3=    atan2((2*cot(beta_3)*(Mn_3^2-1))/(M_3^2*(gam+cos(2*beta_3))+2),1);
T_4=        TR_3*T_3;
gam=        Gamma(T_4);
Mn_4=       sqrt(((gam-1)*Mn_1^2+2)/(2*gam*Mn_1^2-gam+1));
M_4=        Mn_4/sin(beta_3-theta_3);
P_4=        PR_3*P_3;

%fourth shock (internal) 
%theta_4=    theta_1+theta_2+theta_3;
beta_4=     Oblique(M_4,theta_4,T_4);
Mn_4=       M_4*sin(beta_4);
dR_4=       (gam+1)*Mn_4^2/((gam-1)*Mn_4^2+2);
PR_4=       1+2*gam/(gam+1)*(Mn_4^2-1);
TR_4=       PR_4/dR_4;
M_e=        Mn_4/sin(beta_4-theta_4);
P_e=        PR_4*P_4;
T_e=        TR_4*T_4;