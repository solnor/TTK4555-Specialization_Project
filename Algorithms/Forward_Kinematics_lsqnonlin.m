clc; clear all;
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
l_n = [1.73; 3.07; 2.82; 1.29];
l_n = [3.21733641392448;
1.69406945891081;
1.20652780472813;
2.93791763627266];
l_n = [3.21733641392448;
1.69406945891081;
1.20652780472813;
2.93791763627266];
l_n = [3.15444328229476
1.66891906199258
1.20435929655499
2.91702079633575];
n = size(y_n,1);
m = size(l_n, 1);


y_0 = getInitialPoseGuess(a_n,b_n, l_n);
opt = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'OptimalityTolerance',1e-10);
y = lsqnonlin(@(y_i)getPhi(m, a_n, b_n, y_i, l_n).^2, y_0, [], []);
y

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

function phi = getPhi(m, a, b, y, l)
    phi = zeros(m,1);
    for i = 1:m
        phi(i) = generateNu(a(:,i), b(:,i), y(:), l(i));
    end
end

function nu = generateNu(a, b, y, l)
    nu = ((a(1) - y(1) - b(1)*cos(y(3)) + b(2)*sin(y(3)))^2 + (y(2) - a(2) + b(2)*cos(y(3)) + b(2)*sin(y(3)))^2 - l^2)^2;
end