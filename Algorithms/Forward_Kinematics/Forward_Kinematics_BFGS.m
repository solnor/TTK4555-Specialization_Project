clear all; clc;

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
y_i = initialPoseEstimate(a,b,l)

n = size(y_i,1);
m = size(l, 1);

epsilon = 1e-17;

J = Jacobian(a, b, y_i, l);
size(J)
H_0 = eye(3);
phi = getPhi(a, b, y_i, l, 2000, m);
size(phi.')
size(J.')
p_0 = -H_0*J.'*phi
size(p_0)

alpha_bar = 20;
rho = 0.90;

c = 1e-10;

alpha = alpha_bar;
Phi = 0.5*(phi.')*phi;
H_k = eye(3);

i = 0;
BFGS_stop = false;
i_max = 2000;
while ~BFGS_stop && i < i_max

    J = Jacobian(a, b, y_i, l);
    phi = getPhi(a, b, y_i, l, 2000, m);
    F_grad = J.'*phi;
    if norm(F_grad,2) < epsilon
        BFGS_stop = true;
        disp('Stopping BFGS: Condition met.')
        break;
    end
    p_k = -H_k*F_grad;
%     p_k = p_k/norm(p_k,2);

    % Backtracking line search
    Phi = 0.5*(phi.')*phi;
    BLS_stop = false;
    alpha = alpha_bar;
    while ~BLS_stop
        phi_xk_plus_step = getPhi(a, b, y_i+alpha*p_k, l, 2000, m);
        Phi_xk_plus_step = 0.5*(phi_xk_plus_step.')*phi_xk_plus_step;
        if Phi_xk_plus_step <= Phi + c*alpha*(F_grad.')*p_k
            BLS_stop = true;
            break;
        end
        alpha = rho*alpha;
    end
   
    alpha_k = alpha;
    y_inext = y_i + alpha_k*p_k;
    s_k = y_inext - y_i;
    J_knext = Jacobian(a, b, y_inext, l);
    phi_inext = getPhi(a, b, y_inext, l, 2000, m);
    F_grad_inext = J.'*phi_inext;
    gamma_k = F_grad_inext - F_grad;

    H_knext = (eye(3) - rho*s_k*gamma_k.')* ...
              H_k* ...
              (eye(3)-rho*gamma_k*s_k.') + ...
              rho*s_k*(s_k.');
    y_i = y_inext;
    i = i+1;
    if i >= i_max
        disp('Stopping BFGS: Maximum iterations reached')
    end
end
y_i
i