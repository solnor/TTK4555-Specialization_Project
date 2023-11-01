clear all; clc;

% function f_bar getWrench()
a = [-2.0 2.0; 
      2.0 2.0; 
      2.0 0.0; 
     -2.0 0.0].';
% a_n(dimension, cable)
b = [-0.05 0.1;
      0.05 0.1;
      0.05 0.0;
     -0.05 0.0].';
% b_n(dimension, cable)
r = [0;0.7];
x = [r; 0];
% l = [2.61; 1.71; 1.76; 2.65];

m = 4;
n = 3;

l   = inverseKinematics(a, b, x, m)
u   = cableUnitVectors(a, b, x, l, m);
A_T = structureMatrix(u, x(3), b, m);
A_T = A_T(2:4,:)
m_p = 0;
w = [0;-9.81*m_p;0];

f_min = 10*ones(m,1); %[N]
f_max = 200*ones(m,1); %[N]

f_ref = 0.5*(f_min+f_max)
%%
opt = optimoptions('quadprog','Display','none', 'Algorithm', 'interior-point-convex');
f_opt = quadprog(eye(m),-f_ref,[],[],A_T, -w,f_min,f_max,[], opt)
A_T*f_opt
% f_opt = quadprog(eye(m),[],[],[],A_T, -w,f_min,f_max)

