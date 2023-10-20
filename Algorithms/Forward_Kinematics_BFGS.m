clear all; clc;

a_n = [-2.0 2.0; 
        2.0 2.0; 
        2.0 0.0; 
       -2.0 0.0].';
% a_n(dimension, cable)
b_n = [-0.05 0.1;
        0.05 0.1;
        0.05 0.0;
       -0.05 0.0].';
% b_n(dimension, cable)
r_n = [0.5;1];
y_n = [r_n; 0];

n = size(y_n,1);
m = size(l_n, 1);
l_n = inverseKinematics(a_n, b_n, [-0.8; 1.1; -1.0], m);


r_low_i  = zeros(2, m);
r_high_i = zeros(2, m);

for i=1:m
    r_low_i(:, i)  = ...
        a_n(:,i) - (l_n(i) + norm(b_n(:,i), 2))*[1;1];
    r_high_i(:, i) = ...
        a_n(:,i) + (l_n(i) + norm(b_n(:,i), 2))*[1;1];
end

r_low = [max(r_low_i(1,:)); 
         max(r_low_i(2,:))];

r_high = [min(r_high_i(1,:));
          min(r_high_i(2,:));];

r_0 = [0.5*(r_low(1) + r_high(1)); 
       0.5*(r_low(2) + r_high(2)); 
       0.0];
y_k = r_0
% y_k = [0.325; -0.8; 0.6];





epsilon = 1e-7;

J = Jacobian(a_n, b_n, y_k, l_n);
size(J)
H_0 = eye(3);
phi = zeros(4,1);
phi = getPhi(phi, a_n, b_n, y_k, l_n);
size(phi.')
size(J.')
p_0 = -H_0*J.'*phi
size(p_0)


% With H
% p_0 =
% 
%    1.0e+03 *
% 
%    -8.2742
%     8.4856
%    -0.4827

%Without H
% p_0 =
% 
%    1.0e+03 *
% 
%     8.2742
%    -8.4856
%     0.4827

alpha_bar = 1;
rho = 0.8;
% c = 0.25;
c = 1e-3;
% c = 0.1;
% c = 0.1;
alpha = alpha_bar;
Phi = 0.5*(phi.')*phi;
H_k = eye(3);
G = diag([1; 1; 1]);


k = 0;
BFGS_stop = false;
k_max = 1500;
% y_k = [0.9; 0.6; 0];
while ~BFGS_stop && k < k_max
    J = Jacobian(a_n, b_n, y_k, l_n);
    phi = getPhi(phi, a_n, b_n, y_k, l_n);
    F_grad = G.'*J.'*phi;
    if norm(F_grad,2) < epsilon
        BFGS_stop = true;
        disp('Stopping BFGS: Condition met.')
        break;
    end
    p_k = -H_k*F_grad;
%     H_k
%     J
%     F_grad
    % Backtracking line search
    Phi = 0.5*(phi.')*phi;
    BLS_stop = false;
%     alpha_k = 0;
% phi =
% 
%     5.0389
%    51.1053
%     3.4040
%     6.8260
% p_k
    while ~BLS_stop
        phi_xk_plus_step = getPhi(phi, a_n, b_n, y_k+alpha*p_k, l_n);
%         break;
        Phi_xk_plus_step = 0.5*(phi_xk_plus_step.')*phi_xk_plus_step;
        if Phi_xk_plus_step <= Phi + c*alpha*(F_grad.')*p_k
            BLS_stop = true;
%             disp('Found suitable alpha:')
%             disp(alpha)
            break;
        end
        alpha = rho*alpha;
    end
    

    alpha_k = alpha*1.1;
    y_knext = y_k + alpha_k*p_k;
%     y_knext(1:2) = [0.9;0.6];
    s_k = y_knext - y_k;
    J_knext = Jacobian(a_n, b_n, y_knext, l_n);
    phi_knext = getPhi(phi, a_n, b_n, y_knext, l_n);
    F_grad_knext = J.'*phi;
    gamma_k = F_grad_knext - F_grad;

    H_knext = (eye(3) - rho*s_k*gamma_k.')* ...
              H_k* ...
              (eye(3)-rho*gamma_k*s_k.') + ...
              rho*s_k*(s_k.');
    y_k = y_knext
    k = k+1;
%     break;
end
y_k
k

function J = generateNumericJacobian(a_n, b_n, y_k, l_n)
    J_row1 = generateNumericJacobianRow2(a_n(:,1), b_n(:,1), y_k, l_n(1));
    J_row2 = generateNumericJacobianRow2(a_n(:,2), b_n(:,2), y_k, l_n(2));
    J_row3 = generateNumericJacobianRow2(a_n(:,3), b_n(:,3), y_k, l_n(3));
    J_row4 = generateNumericJacobianRow2(a_n(:,4), b_n(:,4), y_k, l_n(4));
    J = [J_row1; J_row2; J_row3; J_row4];
end

function J = Jacobian(a, b, y, l)
    n = size(y,1);
    m = size(a,2);

    J = zeros(m,n);
    for i=1:m
        J(i,:) = JacobianRow(a(:,i), b(:,i), y, l(i));
    end
end

function J_row = JacobianRow(a,b,y_i,l)
    ax = a(1);
    ay = a(2);
    bx = b(1);
    by = b(2);
    x = y_i(1);
    y = y_i(2);
    theta = y_i(3);
    k_spring = 1;
    J_x = ...
    -2*k_spring*(    (ax - x - bx*cos(theta) + by*sin(theta))^2        ...
                 +   (y - ay + by*cos(theta) + bx*sin(theta))^2 - l^2) ...
                 * 2*(ax - x - bx*cos(theta) + by*sin(theta));
 
    J_y = ...
     2*k_spring*(   (ax - x - bx*cos(theta) + by*sin(theta))^2         ...
                 +  (y - ay + by*cos(theta) + bx*sin(theta))^2 - l^2)  ...
                 *2*(y - ay + by*cos(theta) + bx*sin(theta));
 
    J_theta = ...
     2*k_spring*(  2*(by*cos(theta) + bx*sin(theta))                   ...
                    *(ax - x - bx*cos(theta) + by*sin(theta))          ...
                 + 2*(bx*cos(theta) - by*sin(theta))                   ...
                    *(y - ay + by*cos(theta) + bx*sin(theta)))         ...
                   *((ax - x - bx*cos(theta) + by*sin(theta))^2        ...
                   + (y - ay + by*cos(theta) + bx*sin(theta))^2 - l^2);
    J_row = [J_x J_y J_theta];
end


function J_row = generateNumericJacobianRow2(a, b, y, l)
    J_x = -2*((a(1) - y(1) - b(1)*cos(y(3)) + b(2)*sin(y(3)))^2 ...
          + (y(2) - a(2) + b(2)*cos(y(3)) + b(1)*sin(y(3)))^2   ...
          - l^2) ...
          * 2*(a(1) - y(1) - b(1)*cos(y(3)) + b(2)*sin(y(3)));

    J_y = 2*( ( a(1) - y(1) - b(1)*cos(y(3)) + b(2)*sin(y(3)) )^2 ...
            + ( y(2) - a(2) + b(2)*cos(y(3)) + b(1)*sin(y(3)) )^2 ... 
            - l^2)                                                ...
            *2*(y(2) - a(2) + b(2)*cos(y(3)) + b(1)*sin(y(3)));
    J_theta = 2*(2*(b(2)*cos(y(3)) + b(1)*sin(y(3)))*(a(1) - y(1) - b(1)*cos(y(3)) + b(2)*sin(y(3))) + 2*(b(1)*cos(y(3)) - b(2)*sin(y(3)))*(y(2) - a(2) + b(2)*cos(y(3)) + b(1)*sin(y(3))))*((a(1) - y(1) - b(1)*cos(y(3)) + b(2)*sin(y(3)))^2 + (y(2) - a(2) + b(2)*cos(y(3)) + b(1)*sin(y(3)))^2 - l^2);
    J_row = [J_x J_y J_theta];
end


function nu = generateNu2(a, b, y, l)
    nu = ((a(1) - y(1) - b(1)*cos(y(3)) + b(2)*sin(y(3)))^2 + (y(2) - a(2) + b(2)*cos(y(3)) + b(2)*sin(y(3)))^2 - l^2)^2;
end

function phi = getPhi(phi, a, b, y, l)
    for i = 1:size(phi,1)
        phi(i) = generateNu2(a(:,i), b(:,i), y(:), l(i));
    end
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

function Rot = R_z(theta)
    Rot = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
end