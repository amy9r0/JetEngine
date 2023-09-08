close all; clear variables; clc;
R=              287;                            % J/kg/K
g=              9.81;                           % m/s^2
h_PR=           50.01E6;                        % J/kg, methane
rho_v=          100;                            % kg/m^3
f_st=           (36+3*4)/103/(4+4);             % stoichiometric, methane
n=              2;                              % turning load factor
L_burner=       1.2;                            % m
L_inlet=        4;                              % m
theta_1=        5*pi/180;                       % first inlet wedge, optimized for cruise
theta_2=        0.098244067438482;
theta_3=        0.111004110664325;
theta_4=        0.296514640702524;
angle=          10*pi/180;                      % nozzle
tau_b=          1.1;                            % same for ramjet and scramjet
tau_b_tjab=     3;
tau_c=          25;                             % turbojet stage
alt=            35E3;                           % m
% atmosphere calculations
[P_inf,T_inf,rho_inf]=  atmosphere(alt);        % Pa, K
gam_inf=                Gamma(T_inf);
% acceleration stage computation
accel=          20;                             % m/s^2, initial acceleration
v_cruise=       6*sqrt(gam_inf*R*T_inf);        % m/s
t_accel=        v_cruise/accel;                 % s, time to accelerate to cruise
t=              linspace(0,t_accel,61);         % s, time for acceleration
v=              accel*t;                        % m/s, velocity during acceleration
R_accel=        0.5*accel*t(end)^2;             % m, range of acceleration
M_range=        v/sqrt(gam_inf*R*T_inf);
% range calculations
r_turn=         v_cruise^2/g/sqrt(n^2-1);       % m
R_turn=         pi*r_turn;                      % m, range of turn
range=          2500E3*2+R_turn;                % m, total range
R_cruise=       range-R_accel;                  % m, range at cruise
t_cruise=       R_cruise/v_cruise;              % s, duration of cruise
% inlet computation
% three external shocks with constant pressure ratios
% one internal shock
for i=1:length(t)
M=      M_range(i);
if M<3
    burner= 'turbojet';
elseif M>=3 && M<5
    burner= 'ramjet';
else
    burner= 'scramjet';
end
% burner computation, output v_exit/v_inlet, P_exit
switch burner
    case 'turbojet'
        [P_3,T_3,M_3,b_1,b_2,b_3,b_4,P_2i,P_3i,P_4i]=...
            Three_ext_1int(theta_1,theta_2,theta_3,theta_4,3,P_inf,T_inf);
        % geometry/lift/drag of inlet
        [L_1,L_2,L_3,vol_ib,y_i]=...
            Coordinates(theta_1,theta_2,theta_3,b_1,b_2,b_3,b_4,L_inlet,L_burner,0,angle,3);
        P_2i=   P_inf;
        P_3i=   P_inf;
        P_4i=   P_inf;
        T_0_3=  T_3*(1+(Gamma(T_3)-1)/2*M_3^2);
        [ve_v0(i),P_5]= Turbojet(T_inf,P_inf,h_PR,f_st,tau_b_tjab,M,tau_c);
    case 'ramjet'
        [P_3,T_3,M_3,b_1,b_2,b_3,b_4,P_2i,P_3i,P_4i]=...
            Three_ext_1int(theta_1,theta_2,theta_3,theta_4,M,P_inf,T_inf);
        % geometry/lift/drag of inlet
        [L_1,L_2,L_3,vol_ib,y_i]=...
            Coordinates(theta_1,theta_2,theta_3,b_1,b_2,b_3,b_4,L_inlet,L_burner,0,angle,M);
        T_0_3=  T_3*(1+(Gamma(T_3)-1)/2*M_3^2);
        % setup IC for Shapiro's equations and Runge-Kutta 4th Order method
        IC= [T_0_3;M_3^2;T_3;P_3];
        theta=  45;
        [T_0_4,M_sq_4,T_4(i),P_4]= Shapiro(0,0.05,IC,tau_b,L_burner,theta);
        M_4=    sqrt(M_sq_4);
        if isnan(T_0_4) || isnan(M_sq_4) || isnan(T_4(i)) || isnan(P_4)
            disp('Burner has reached NaN values.')
        end
    case 'scramjet'
        [P_3,T_3,M_3,b_1,b_2,b_3,b_4,P_2i,P_3i,P_4i]=...
            Three_ext_1int(theta_1,theta_2,theta_3,theta_4,M,P_inf,T_inf);
        % geometry/lift/drag of inlet
        [L_1,L_2,L_3,vol_ib,y_i]=...
            Coordinates(theta_1,theta_2,theta_3,b_1,b_2,b_3,b_4,L_inlet,L_burner,1,angle,M);
        T_0_3=  T_3*(1+(Gamma(T_3)-1)/2*M_3^2);
        % setup IC for Shapiro's equations and Runge-Kutta 4th Order method
        IC= [T_0_3;M_3^2;T_3;P_3];
        theta=  5;
        [T_0_4,M_sq_4,T_4(i),P_4]= Shapiro(0,0.05,IC,tau_b,L_burner,theta);
        if isnan(T_0_4) || isnan(M_sq_4) || isnan(T_4(i)) || isnan(P_4)
            disp('Burner has reached NaN values.')
        end
        M_4=    sqrt(M_sq_4);
        if M_3<=1 || M_4<=1
            disp('The scramjet has reached subsonic speed.')
        end
end
switch burner
    case {'scramjet','ramjet'}
        % Prandtl-Meyer expansion wave
        nu=             angle+sqrt((Gamma(T_4(i))+1)/(Gamma(T_4(i))-1))*...
                        atan2(sqrt((Gamma(T_4(i))-1)/(Gamma(T_4(i))+1)*(M_4^2-1)),1)...
                        -atan2(sqrt(M_4^2-1),1);
        f=              @(M)(sqrt((Gamma(T_4(i))+1)/(Gamma(T_4(i))-1))*...
                        atan2(sqrt((Gamma(T_4(i))-1)/(Gamma(T_4(i))+1)*(M^2-1)),1)...
                        -atan2(sqrt(M^2-1),1)-nu);
        [M_5,error(i)]= FalsePosition(1,7,1E-7,f,1E5);
        % compute exit properties
        T_5=            T_0_4/(1+(Gamma(T_4(i))-1)/2*M_5^2);
        gam=            Gamma(T_5);
        P_5=            P_4*(1+(Gamma(T_4(i))-1)/2*M_4^2)^(Gamma(T_4(i))/(Gamma(T_4(i))-1))...
                        /(1+(gam-1)/2*M_5^2)^(gam/(gam-1));
        ve_v0(i)=       M_5/M*sqrt(gam*T_5/gam_inf/T_inf);
end
if ve_v0<=1
   disp('The engine is unstarted.') 
end
L_i=                        L_1*P_2i*cos(theta_1)+L_2*P_3i*cos(theta_1+theta_2)+...
                            L_3*P_4i*cos(theta_1+theta_2+theta_3);
D_i=                        L_1*P_2i*sin(theta_1)+L_2*P_3i*sin(theta_1+theta_2)+...
                            L_3*P_4i*sin(theta_1+theta_2+theta_3);
% mass flow rate, total initial area assumed to catch air mass
m0(i)=      -y_i*rho_inf*M*sqrt(gam_inf*R*T_inf);               % kg/s
% characteristics
F_m0=       ((1+f_st)*ve_v0(i)-1)*M*sqrt(gam_inf*R*T_inf);      % N*s/kg
I_sp(i)=    F_m0/g/f_st;                                        % s
L_nozzle=   -y_i*tan(pi/2-angle);                               % m
vol=        vol_ib+0.5*-y_i*L_nozzle;                           % m^3
weight=     vol*rho_v;                                          % kg
T_surf=     sqrt(L_nozzle^2+y_i^2);                             % m
L_n=        P_5*T_surf*sin(angle);                              % lift per width of nozzle
Thrust(i)=  P_5*T_surf*cos(angle)+F_m0*m0(i)-D_i;               % thrust per width
phi_e(i)=   D_i/(F_m0*m0(i)+P_5*T_surf*cos(angle));             % installation loss
L(i)=       L_n+L_i;                                            % lift per width
L_D(i)=     L(i)/D_i;                                           % lift to drag
L_W(i)=     L(i)/weight;                                        % lift to weight
end
m_f_accel=  CompositeSimpson(0,t_accel,f_st*m0,length(t)-1);    % kg
m_f_cruise= f_st*m0(end)*t_cruise;                              % kg
m_p=        Thrust(end)/accel-weight;                           % kg
m_f=        m_f_accel+m_f_cruise;                               % kg
m_i=        weight+m_f+m_p;                                     % initial mass of aircraft
pi_f=       m_f/m_i;                                            % fuel mass fraction
figure
plot(t,ve_v0,'.k'), title('Velocity Ratio over Time where Acceleration is 20 m/s^2.') 
xlabel('t (s)'), ylabel('v_e/v_0 (Non-Dimensional)'), grid
figure
plot(t,I_sp,'.k'), title('Specific Impulse over Time where Acceleration is 20 m/s^2.')
xlabel('t (s)'), ylabel('I_{sp} (s)'), grid
figure
plot(t,Thrust/1E3,'.k'), title('Thrust per Width over Time where Acceleration is 20 m/s^2.')
xlabel('t (s)'), ylabel('T/w (kN/m)'), grid
figure
plot(t,phi_e,'.k'), title('Installation Loss over Time where Acceleration is 20 m/s^2.')
xlabel('t (s)'), ylabel('\phi_e (Non-Dimensional)'), grid
figure
plot(t,L_D,'.k'), title('Lift to Drag over Time where Acceleration is 20 m/s^2.')
xlabel('t (s)'), ylabel('L/D (Non-Dimensional)'), grid
figure
plot(t,L_W,'.k'), title('Lift to Weight over Time where Acceleration is 20 m/s^2.')
xlabel('t (s)'), ylabel('L/W (Non-Dimensional)'), grid
figure
plot(t,L/1E3,'.k'), title('Lift over Time where Acceleration is 20 m/s^2.')
xlabel('t (s)'), ylabel('L (kN)'), grid
figure
plot(t,f_st*m0,'.k'), title('Mass Flow-rate of Fuel per Width over Time where Acceleration is 20 m/s^2.')
xlabel('t (s)'), ylabel('dm_f/dt/w (kg/s/m)'),grid
figure
plot(t,T_4,'.k'), title('Exit temperature of Burner over Time where Acceleration is 20 m/s^2.')
xlabel('t (s)'), ylabel('T_4 K'), grid