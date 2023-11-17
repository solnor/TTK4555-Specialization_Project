clear all; clc;

% function f_bar getWrench()
a = [-2.0 2.0; 
      2.0 2.0; 
      2.0 0.0; 
     -2.0 0.0].';
% a_n(dimension, cable)
b = [-0.05  0.05;
      0.05  0.05;
      0.05 -0.05;
     -0.05 -0.05].';
% b_n(dimension, cable)

r = [0;1];

x = [r; 0];
x_0 = [x(1); 0; x(2); 0; x(3);0];
x_d = [0.0;0.5;0];
x_ss_0 = [x(1); 0; x(2); 0; x(3);0];

m = 4;
n = 3;

l_0   = inverseKinematics(a, b, x, m)
u   = cableUnitVectors(a, b, x, l_0, m);
A_T = structureMatrix(u, x(3), b, m);
A_T = A_T(2:4,:)
m_p = 0.1;
w = [0;-9.81*m_p;0];

f_min = 2*ones(m,1); %[N]
f_max = 60*ones(m,1); %[N]

f_ref = 0.5*(f_min+f_max)

f_0 = 15*ones(4,1);

k_spring = diag(200*ones(m,1))

k_spring_inv = inv(k_spring);

A_p = [0  1   0  0   0  0;
       0 -2.0 0  0   0  0;
       0  0   0  1   0  0;
       0  0   0 -2.0 0  0;
       0  0   0  0   0  1;
       0  0   0  0   0 -2.0;];
A_p = [0  1   0  0   0  0;
       0 -2.0 0  0   0  0;
       0  0   0  1   0  0;
       0  0   0 -40.0 0  0;
       0  0   0  0   0  1;
       0  0   0  0   0 -20.0;];

I_zz = 2*(0.05*50e-3+0.05*50e-3)*100;

B_p = [0     0     0;
       1/m_p 0     0;
       0     0     0;
       0     1/m_p 0;
       0     0     0;
       0     0     1/I_zz;];

B_p_d = [-0.2  0    0;
          0    0    0;
          0   -0.2  0;
          0    0    0;
          0    0   -0.2;
          0    0    0];


r_winch = 37.5e-3;
r_winch = 60.0e-3;
psi = l_0/r_winch
l_0(1)
f_0(1)
k_spring(1,1)
l_0(1)-f_0(1)*k_spring_inv(1,1)
l_0(1)-f_0(1)*k_spring_inv(1,1)
theta_r_0 = (l_0(1)-f_0(1)*k_spring_inv(1,1))/r_winch;
l_0(1)-f_ref(1)*k_spring_inv(1,1)
theta_r_0 = (l_0(1)-f_ref(1)*k_spring_inv(1,1))/r_winch;


% MarLev position estimation constants
tau       = 1e-6;
epsilon_1 = 1e-17;
epsilon_2 = 1e-17;
xi        = 2;

motor_init;