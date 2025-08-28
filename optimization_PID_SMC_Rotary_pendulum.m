syms theta alpha dtheta dalpha ddtheta ddalpha  real
syms mp Lr Lp Jr Jp g Dr Dp tau real
syms t phi
c1=4.6210
c2=4.1349


mp = 0.024;        % Mass of pendulum (kg)
Lp = 0.129;       % Length of pendulum (m)
Lr = 0.085;       % Length of rotary arm (m)
Jp = 3.3282*10^-5;       % Pendulum inertia (kg.m^2)
Jr = 5.71*10^-5;        % Rotary arm inertia (kg.m^2)
Dp = 0.0005;       % Pendulum damping
Dr = 0.0015;        % Rotary arm damping
g = 9.8;



% M(q)
M11 = Jr + mp*Lr^2 + (1/4)*mp*Lp^2*sin(alpha)^2;
M12 = (1/2)*mp*Lp*Lr*cos(alpha);
M = [M11, M12; M12, Jp + (1/4)*mp*Lp^2];
Minv=inv(M)
% C(q)
C1 = -0.5*mp*Lp^2*sin(2*alpha)*dalpha^2 - mp*Lp*Lr*sin(alpha)*dalpha*dtheta;
C2 = 0.5*mp*Lp*Lr*sin(alpha)*dtheta^2;
C = [C1; C2];
% G(q)
G = [0; -0.5*mp*g*Lp*sin(alpha)];
% Damping
D = [Dr*dtheta; Dp*dalpha];



 

K1 = Jp + (1/4) * mp * Lp^2;
K2 = (1/2) * mp * Lp * Lr;
K3 = Jr + mp * Lr^2;
JT = K1 *  K3 - K2^2;

% State-Space Matrices
A = (1/JT) * [
    0, 0, JT, 0;
    0, 0, 0, JT;
    0,(1/4)*mp^2*Lp^2*Lr*g, -K1*Dr, K2*Dp;
   0, -(1/2)*mp*Lp*g*K3, (1/2)*mp*Lp*Lr*Dr, -K3*Dp
];

B = (1/JT) * [
    0;
    0;
    K1;
    -K2
];
BB=B'

C = [1 0 0 0;
     0 1 0 0];

D = [0; 0];

 
 eq1 = (mp*Lr^2 + Jr + 0.25*mp*Lp^2 - 0.25*mp*Lp^2*cos(alpha)^2)*ddtheta ...
    - 0.5*mp*Lp*Lr*cos(alpha)*ddalpha ...
    + 0.5*mp*Lp^2*sin(alpha)*cos(alpha)*dtheta*dalpha ...
    + 0.5*mp*Lp*Lr*sin(alpha)*dalpha^2 ...
    + Dr*dtheta - tau == 0;
 
eq2 = 0.5*mp*Lp*Lr*cos(alpha)*ddtheta ...
    + (Jp + 0.25*mp*Lp^2)*ddalpha ...
    - 0.25*mp*Lp^2*sin(alpha)*cos(alpha)*dtheta^2 ...
    + 0.5*mp*Lp*g*sin(alpha) ...
    + Dp*dalpha == 0;
 
sol = solve([eq1, eq2], [ddtheta, ddalpha])
 
 
x0 = [0; 5*pi/180; 0; 0];



t=10
pha1=0.1*sin(pi*t);
pha2=0.1*sin(pi/2*t);
kp=15.94461;
ki=1.9406;
kd=1.6452;

 % V = phi_d_ddot + obj.c1 * e + obj.c2 * e_dot;
   %u_smc = (V - g_func(phi, phi_dot)) / h_func(phi, phi_dot);
%slid mode contlrol parmater 

 GSMC=-Minv *[C+G];
 HSMC=Minv*[0 ;1]
 
