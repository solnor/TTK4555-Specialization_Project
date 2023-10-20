clear all; clc;
% a(dimension, cable)
a = [-2.0 2.0; 
      2.0 2.0; 
      2.0 0.0; 
     -2.0 0.0].';

b = [-0.05 0.1;
      0.05 0.1;
      0.05 0.0;
     -0.05 0.0].';
% b_n(dimension, cable)


l   = inverseKinematics(a, b, [-0.8; 1.1; -1.0], size(a, 2));
y_i = initialPoseEstimate(a,b,l);

n = size(y_i, 1);
m = size(a, 2);

J = Jacobian(a, b, y_i, l, 1);
A = J.'*J;
phi = getPhi(a, b, y_i, l, m);
g = J.'*phi;


tau       = 1e-6;
mu        = tau*max(diag(A));
epsilon_1 = 1e-17;
epsilon_2 = 1e-17;
xi        = 2;

iter = 0;
max_iter = 100;

stop = false;
while ~stop && iter < max_iter
    iter = iter+1;

    h = (A + mu*eye(3))\(-g);
    
    if norm(h,2) <= epsilon_2*(norm(y_i,2)+epsilon_2)
        stop = true;
        disp(epsilon_2*(norm(y_i,2)+epsilon_2));
        disp('Stopping algorithm')
    else
        y_new = y_i + h;

        % Gain ratio
        Phi = 0.5*(phi.')*phi;
        phi_x_plus_h = getPhi(a, b, y_new, l, m);
        
        Phi_x_plus_h = 0.5*(phi_x_plus_h.')*phi_x_plus_h;
        rho = (Phi - Phi_x_plus_h)/(0.5*(h.')*(mu*h - g));

        if rho > 0
            y_i = y_new;
            J = Jacobian(a, b, y_i, l, 1);
            A = J.'*J;
            phi = getPhi(a, b, y_i, l, m);
            g = J.'*phi;
            cond = norm(g,2) < epsilon_1;
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
    end
end
cond = vpa(norm(g,2)) < epsilon_1;
iter
y_i



