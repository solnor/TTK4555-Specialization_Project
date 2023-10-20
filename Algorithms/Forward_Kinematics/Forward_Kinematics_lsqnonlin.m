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


% y_0 = getInitialPoseGuess(a,b, l);
opt = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'OptimalityTolerance',1e-10);
y = lsqnonlin(@(y_i)getPhi(a, b, y_i, l, 2000, m).^2, y_0, [], [], opt);
y