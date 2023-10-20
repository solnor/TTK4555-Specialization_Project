clc; clear all;
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


l = inverseKinematics(a, b, [-0.8; 1.1; -1.0], size(a, 2));
y_0 = initialPoseEstimate(a, b, l)

n = size(y_0,1);
m = size(l, 1);

opt = optimoptions('fmincon',            ...
                   'Algorithm',   'interior-point', ...
                   'MaxFunEvals', 2000);
[z, ~, ~, output]   = fmincon(@(y_i)objfunc(y_i, a, b, l), y_0, ...
                              [],   [], [], [],[],[],[],opt);
z
output

function f = objfunc(y_i, a_n, b_n, l_n)
    x = y_i(1);
    y = y_i(2);
    theta = y_i(3);
    n1 = (norm(a_n(:,1) - [x;y] - R_z(theta)*b_n(:,1), 2)^2 - l_n(1)^2)^2;
    n2 = (norm(a_n(:,2) - [x;y] - R_z(theta)*b_n(:,2), 2)^2 - l_n(2)^2)^2;
    n3 = (norm(a_n(:,3) - [x;y] - R_z(theta)*b_n(:,3), 2)^2 - l_n(3)^2)^2;
    n4 = (norm(a_n(:,4) - [x;y] - R_z(theta)*b_n(:,4), 2)^2 - l_n(4)^2)^2;

    f = 0.5*(n1+n2+n3+n4)^2;
end

function R = R_z(theta)
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
end
