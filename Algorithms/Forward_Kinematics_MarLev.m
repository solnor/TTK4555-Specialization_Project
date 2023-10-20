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
l_n = [2.61; 1.71; 1.76; 2.65];
l_n = [1.67; 2.34; 2.78; 2.19];
l_n = [1.63; 2.36; 2.78; 2.17];
l_n = [1.73; 3.07; 2.82; 1.29];
l_n = [3.21733641392448;
1.69406945891081;
1.20652780472813;
2.93791763627266];
l_n = [3.13249102153542;
1.67107749670684;
1.20933866224478;
2.91247317584214];
l_n = [3.14565683386651;
1.66870828876108;
1.20579216065453;
2.91493247924206];
l_n = [3.15444328229476;
1.66891906199258;
1.20435929655499;
2.91702079633575;];
l_n = [1.67377392891419;
2.85811042163665;
2.76637985913162;
1.48926443864987;];

n = size(y_n,1);
m = size(l_n, 1);


r_low_i  = zeros(2, m);
r_high_i = zeros(2, m);




% for i=1:m
%     r_low_i(:, i)  = ...
%         a_n(:,i) - (l_n(i) + norm(b_n(:,i), 2))*[1;1];
%     r_high_i(:, i) = ...
%         a_n(:,i) + (l_n(i) + norm(b_n(:,i), 2))*[1;1];
% end
% 
% r_low = [max(r_low_i(1,:)); 
%          max(r_low_i(2,:))];
% 
% r_high = [min(r_high_i(1,:));
%           min(r_high_i(2,:));];
% 
% r_0 = [0.5*(r_low(1) + r_high(1)); 
%        0.5*(r_low(2) + r_high(2)); 
%        0.0];
% y_i = r_0


% y_i = [0.325; -0.8; 0.6];



% 
% phi = zeros(4,1);
% for i = 1:size(phi,1)
%     phi(i) = generateNu2(a_n(:,i), b_n(:,i), y_i(:), l_n(i));
% end
% g = J.'*phi;

% h = (A+mu*eye(3))\(-g);
% % h = [0.01; 0.01; 0.01];
% Phi = 0.5*(phi.')*phi
% 
% phi_x_plus_h = zeros(4,1);
% for i = 1:size(phi_x_plus_h,1)
%     phi_x_plus_h(i) = generateNu2(a_n(:,i), b_n(:,i), y_i(:)+h(:), l_n(i));
% end
% Phi_x_plus_h = 0.5*(phi_x_plus_h.')*phi_x_plus_h
% Phi-Phi_x_plus_h
% denominator = 0.5*(h.')*(mu*h-g)
% rho = (Phi-Phi_x_plus_h)/denominator
% % Phi_taylor = Phi + (h.')*(J.')*phi + 0.5*(h.')*(J.')*J*h

l_n = inverseKinematics(a_n, b_n, [-0.8; 1.1; -1.0], m)
y_i = getInitialPoseGuess(a_n,b_n,l_n)

phi          = zeros(4,1);
phi_x_plus_h = zeros(4,1);

J = generateNumericJacobian(a_n, b_n, y_i, l_n);
A = J.'*J;
phi = getPhi(a_n, b_n, y_i, l_n, m);
g = J.'*phi;

G = diag([1; 1; 1]);


tau       = 1e-6;
mu        = tau*max(diag(A));
epsilon_1 = 10e-17;
epsilon_2 = 10e-17;
xi        = 2;

iter = 0;
max_iter = 1000;

stop = false;
while ~stop && iter < max_iter
    iter = iter+1;

%     J = generateNumericJacobian(a_n, b_n, y_i, l_n);
%     A = J.'*J;
%     phi = getPhi(phi, a_n, b_n, y_i, l_n);
%     g = J.'*phi;
    h = (A + mu*eye(3))\(-g);
    
    if norm(h,2) <= epsilon_2*(norm(y_i,2)+epsilon_2)
        stop = true;
        disp(epsilon_2*(norm(y_i,2)+epsilon_2));
        disp('Stopping algorithm')
    else
        y_new = y_i + h;

        % Gain ratio
%         Phi = 0.5*(phi.')*phi;
        Phi = norm(phi)^2;
        phi_x_plus_h = getPhi(a_n, b_n, y_new, l_n, m);
        
%         Phi_x_plus_h = 0.5*(phi_x_plus_h.')*phi_x_plus_h;
        Phi_x_plus_h = norm(phi_x_plus_h)^2;
        rho = (Phi - Phi_x_plus_h)/(0.5*(h.')*(mu*h - g));

        if rho > 0
            y_i = y_new;
            J = generateNumericJacobian(a_n, b_n, y_i, l_n);
            A = J.'*J;
            phi = getPhi(a_n, b_n, y_i, l_n, m);
            g = J.'*phi;
            cond = norm(g,2) < epsilon_1;
%             cond = max(abs(g)) < epsilon_1;
            if cond
                disp('Stopping algorithm 2')
                stop = cond;
            end
            mu = mu*max([1/3, 1-(2*rho-1)^3]);
            xi = 2;
        else
            mu = mu*xi;
            xi = 2*xi;
        end
        % J = generateNumericJacobian(a_n, b_n, y_i, l_n);
        
%         for i = 1:size(phi,1)
%             phi(i) = generateNu2(a_n(:,i), b_n(:,i), y_i(:), l_n(i));
%         end
        % phi = getPhi(phi, a_n, b_n, y_i, l_n);
        % g = J.'*phi;
        % h = (A + mu*eye(3))\(-g);
        
    end
end
cond = vpa(norm(g,2)) < epsilon_1
iter
y_i
rad2deg(y_i(3))
disp("Done")


function phi = getPhi(a, b, y, l,m)
    phi = zeros(m,1);
    for i = 1:size(phi,1)
        phi(i) = generateNu(a(:,i), b(:,i), y(:), l(i));
    end
end

function phi = getPhi2(phi, a, b, y, l)
    for i = 1:size(phi,1)
        phi(i) = generateNu4(a(:,i), b(:,i), y(:), l(i));
    end
end

function J = generateNumericJacobian(a_n, b_n, y_i, l_n)
%     J_row1 = generateNumericJacobianRow2(a_n(:,1), b_n(:,1), y_i, l_n(1));
%     J_row2 = generateNumericJacobianRow2(a_n(:,2), b_n(:,2), y_i, l_n(2));
%     J_row3 = generateNumericJacobianRow2(a_n(:,3), b_n(:,3), y_i, l_n(3));
%     J_row4 = generateNumericJacobianRow2(a_n(:,4), b_n(:,4), y_i, l_n(4));
%     J_row1 = JacobianRow(a_n(:,1), b_n(:,1), y_i, l_n(1));
%     J_row2 = JacobianRow(a_n(:,2), b_n(:,2), y_i, l_n(2));
%     J_row3 = JacobianRow(a_n(:,3), b_n(:,3), y_i, l_n(3));
%     J_row4 = JacobianRow(a_n(:,4), b_n(:,4), y_i, l_n(4));
    n = size(y_i,1);
    m = size(a_n,2);

    J = zeros(m,n);
    for i=1:m
        J(i,:) = JacobianRow2(a_n(:,i), b_n(:,i), y_i, l_n(i));
    end
%     J_row1 = JacobianRow2(a_n(:,1), b_n(:,1), y_i, l_n(1));
%     J_row2 = JacobianRow2(a_n(:,2), b_n(:,2), y_i, l_n(2));
%     J_row3 = JacobianRow2(a_n(:,3), b_n(:,3), y_i, l_n(3));
%     J_row4 = JacobianRow2(a_n(:,4), b_n(:,4), y_i, l_n(4));
%     J = [J_row1; J_row2; J_row3; J_row4];
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

function J_row = JacobianRow2(a,b,y_i,l)
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


function J_row = JacobianRow(a,b,y_i,l)
% a(1:2)
ax = a(2);
ay = a(1);
bx = b(1);
by = b(2);
x = y_i(1);
y = y_i(2);
theta = y_i(3);
J_x = ...
    2*x - 2*ax + 2*bx*cos(theta) - 2*by*sin(theta);
J_y = ...
    2*y - 2*ay + 2*by*cos(theta) + 2*bx*sin(theta);
J_theta = ...
    2*(by*cos(theta) + bx*sin(theta))*(ax - x - bx*cos(theta) + by*sin(theta)) + 2*(bx*cos(theta) - by*sin(theta))*(y - ay + by*cos(theta) + bx*sin(theta));
J_row = [J_x J_y J_theta];
end

function nu = generateNu(a, b, y, l)
    k_spring = 1;
    nu = k_spring*( ...
            (a(1) - y(1) - b(1)*cos(y(3)) + b(2)*sin(y(3)))^2 ...
          + (y(2) - a(2) + b(2)*cos(y(3)) + b(1)*sin(y(3)))^2 - l^2 ...
         )^2;
end

function nu = generateNu4(a,b,y_i,l)
    ax = a(2);
    ay = a(1);
    bx = b(1);
    by= b(2);
    x = y_i(1);
    y = y_i(2);
    theta = y_i(3);
    nu = (ax - x - bx*cos(theta) + by*sin(theta))^2 + (y - ay + by*cos(theta) + bx*sin(theta))^2 - l^2;
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
    y_0 = [r_0; 0.0];
end

% function nu = generateNu3(a,b,y,l)
%     a
%     y
%     b
%     nu = ((a-y(1:2)-R_z(y(3))*b).'*(a-y(1:2)-R_z(y(3))*b) - l)^2;
% end

function Rot = R_z(theta)
    Rot = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
end