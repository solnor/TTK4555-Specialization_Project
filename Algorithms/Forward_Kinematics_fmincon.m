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
l_n = [3.21733641392448
1.69406945891081
1.20652780472813
2.93791763627266];
n = size(y_n,1);
m = size(l_n, 1);


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
       0];
y_0 = r_0;
opt = optimoptions('fmincon',            ...
                   'Algorithm',   'interior-point', ...
                   'MaxFunEvals', 8000);
[z, ~, ~, output]   = fmincon(@(y_i)objfunc(y_i, a_n, b_n, l_n), y_0, ...
                              [],   [], [], [],[],[],[],opt);
z
output

function f = objfunc(y_i, a_n, b_n, l_n)
    x = y_i(1);
    y = y_i(2);
    theta = y_i(3);
    n1 = (vpa(norm(a_n(:,1) - [x;y] - R_z(theta)*b_n(:,1), 2)^2) - l_n(1)^2)^2;
    n2 = (vpa(norm(a_n(:,2) - [x;y] - R_z(theta)*b_n(:,2), 2)^2) - l_n(2)^2)^2;
    n3 = (vpa(norm(a_n(:,3) - [x;y] - R_z(theta)*b_n(:,3), 2)^2) - l_n(3)^2)^2;
    n4 = (vpa(norm(a_n(:,4) - [x;y] - R_z(theta)*b_n(:,4), 2)^2) - l_n(4)^2)^2;

    f = cast(0.5*(n1+n2+n3+n4)^2, 'double');
end

function Rot = R_z(theta)
    Rot = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
end
