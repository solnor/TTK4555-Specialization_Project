function A_T = calculateStructureMatrix(u, theta, b, m)
    crossProd = zeros(3, m);
    for i=1:m
        crossProd(:,i) = cross([0;R_z(theta)*b(:,i);], [0;u(:,i)]);
    end
    A_T = [zeros(1,m);
           u;
           crossProd];
end

function R = R_z(theta)
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
end