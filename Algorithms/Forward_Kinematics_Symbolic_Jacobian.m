
syms ax ay      real
syms x y theta  real
syms bx by      real
syms l          real
a = [ax; ay];
b = [bx; by];
r = [x; y];
% l = [lx; ly];
% norm(a-r-R_z(theta)*b)^2-l^2
k_spring = 2000;
nu = 2000*(simplify(norm(a-r-R_z(theta)*b,2)^2-l^2))^2

J_11 = jacobian(nu,x)
J_12 = jacobian(nu,y)
J_13 = jacobian(nu,theta)
J_1  = [J_11 J_12 J_13];

function R = R_z(theta)
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
end

% function nu = generateNu4(a,b,y,l)
%     nu = 
% end