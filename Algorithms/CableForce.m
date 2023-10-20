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
r = [0.5;1];
x = [r; 0];
l = [2.61; 1.71; 1.76; 2.65];
l = [3.14565683386651;
1.66870828876108;
1.20579216065453;
2.91493247924206];
l = [3.15444328229476;
1.66891906199258;
1.20435929655499;
2.91702079633575];
l = [3.13249102153542;
1.67107749670684;
1.20933866224478;
2.91247317584214];
l = [1.67107749670684;
3.13249102153542;
2.91247317584214;
1.20933866224478];
l= [2.23886131772381;
2.23886131772381;
2.10772389083580;
2.10772389083580];

x = [0.9;	0.6;	0.0];
% x = [0.900000000000001	0.600000000000001	-0.270000000000000].';

% l_ik - l;



% l_1 = a(:,1) - x(1:2) - R_z(x(3))*b(:,1);
% u_1 = l_1/l_ik(1);

% norm(l_1/l_ik(1),2);

k_spring = 2000;
l_ik = inverseKinematics(a,b,x,4);
f_bar = calculateCableForces(l_ik, l, k_spring);
u = calculateCableUnitVectors(a,b,x,l_ik,l,4);
A_T = calculateStructureMatrix(u, x(3), b, 4);
A_T_MPinv = pinv(A_T)
% size(A_T)
% size(A_T_MPinv)
w_0   = [0; 0; -9.81; 0; 0; 0];
A_T_MPinv*w_0;
null(A_T);

f_min = ones(4,1)*10; %[N]
f_max = ones(4,1)*50; %[N]
f_ref = (f_min + f_max)/2;
f_CF = f_ref - A_T_MPinv*(w_0+A_T*f_ref)
% A_T_MPinv*A_T
f_CF = -A_T_MPinv*w_0 + (eye(4)-A_T_MPinv*A_T)*f_ref
delta_f = f_CF-f_min;
f_CFE = -A_T_MPinv*w_0 + (eye(4)-A_T_MPinv*A_T)*delta_f
% f_v = (A_T*A_T.')\(-w_0-A_T*f_ref)
% F = [1;2;3;4];
% A_T*F;

%%
f_min = ones(4,1)*10; %[N]
f_max = ones(4,1)*40; %[N]
k_spring = 2000;
l_ik = inverseKinematics(a,b,x,4);
% f_bar = calculateCableForces(l_ik, l, k_spring);
u = calculateCableUnitVectors(a,b,x,l_ik,l,4);
A_T = calculateStructureMatrix(u, x(3), b, 4);
A_T = A_T(2:4,:);
A_T_MPinv = pinv(A_T);
h = null(A_T);
w_0   = [0; 0; -9.81; 0; 0; 0];
w_0 = w_0(2:4);
% for i=1:4
%     f_0 = A_T_MPinv*w_0;
%     (f_min(i) + f_0(i))/h(i)
% end
f_0 = A_T_MPinv*w_0;
lambda_l = max((f_min + f_0)./h)
lambda_h = min((f_max + f_0)./h)
sol_exists = lambda_l < lambda_h;
if sol_exists
    disp('Solution exists')
else
    disp('No solution exists')
end
f_ss = f_0 + 0.5*(lambda_h+lambda_l)*h

f_ref = (f_min + f_max)/2;
f_CF = f_ref - A_T_MPinv*(w_0+A_T*f_ref)
delta_f = f_CF-f_min
delta_f = ones(4,1)*max(delta_f)
f_CFE = -A_T_MPinv*w_0 + (eye(4)-A_T_MPinv*A_T)*(delta_f-5*ones(4,1))
ith = 1;
A_T_marked = [A_T(:,1:ith-1) A_T(:,ith+1:end)]
f_marked = [f_CFE(1:ith-1); f_CFE(ith+1:end)]
w_marked = f_min(1)*A_T(:,ith)
A_T_MPinv_marked = pinv(A_T_marked);
f_CFE = -A_T_MPinv_marked*w_marked + (eye(3)-A_T_MPinv_marked*A_T_marked)*f_marked
%%
% phi_OS(a,b,x,l,k_spring, 4)

x_0 = getInitialPoseGuess(a, b, l);
opt = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'OptimalityTolerance',1e-5, 'Display','off');
% opt = optimoptions('lsqnonlin', 'Algorithm','trust-region-reflective', 'OptimalityTolerance',1e-10, 'Display','iter');
[y, ~, ~, ~, output]  = lsqnonlin(@(x)phi_OS(a,b,x,l,k_spring, 4), x_0, [], [], opt);
y
output

function phi = phi_OS(a, b, x, l_meas, k_spring, m)
    
    l_ik  = inverseKinematics(a, b, x, m);
    f_bar = calculateCableForces(l_ik, l_meas, k_spring);
    u     = calculateCableUnitVectors(a, b, x, l_ik, l_meas, m);
    A_T   = calculateStructureMatrix(u, x(3), b, m);
%     W     = diag([2; 1; 1; 1; 1; 1]);
    w_0   = [0; 0; -9.81; 0; 0; 0];
    W = diag(1*ones(6,1));
    phi = W*(A_T*f_bar-w_0).^2;
end

function y_0 = getInitialPoseGuess(a, b, l)
    m = size(a,2);
    r_low_i  = zeros(2, m);
    r_high_i = zeros(2, m);
    
    for i=1:m
        r_low_i(:, i)  = ...
            a(:,i) - (l(i) + norm(b(:,i), 2))*[1;1];
        r_high_i(:, i) = ...
            a(:,i) + (l(i) + norm(b(:,i), 2))*[1;1];
    end
    
    r_low = [max(r_low_i(1,:)); 
             max(r_low_i(2,:))];
    
    r_high = [min(r_high_i(1,:));
              min(r_high_i(2,:));];
    
    r_0 = [0.5*(r_low(1) + r_high(1)); 
           0.5*(r_low(2) + r_high(2))];
    y_0 = [r_0; 0.20];
end

function A_T = calculateStructureMatrix(u, theta, b, m)
    crossProd = zeros(3, m);
    for i=1:m
%         cp_ = cross([R_z(theta)*b(:,i); 0], [u(:,i);0]);
%         crossProd(1,i) = cp_(3,:);
        crossProd(:,i) = cross([0;R_z(theta)*b(:,i);], [0;u(:,i)]);
    end
%     disp("Hey");
%     size(u)
    A_T = [zeros(1,m);
           u;
           crossProd];
end

function u = calculateCableUnitVectors(a, b, x, l_ik, l_meas, m)
    l = a - x(1:2) - R_z(x(3))*b;
    u = zeros(2,4);
    for i=1:m
        u(:,i) = l(:,i)/l_ik(i);
    end
end

function f_bar = calculateCableForces(l_ik, l_meas, k)
    f_bar = k*(l_ik - l_meas);
end

function l = inverseKinematics(a, b, pose, m)
    x     = pose(1);
    y     = pose(2);
    theta = pose(3);

    l = zeros(m,1);
    for i=1:m
        l(i) = norm(a(:,i)-[x;y]-R_z(theta)*b(:,i),2);
    end
end


function R = R_z(theta)
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
end