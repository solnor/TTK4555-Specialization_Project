function u = calculateCableUnitVectors(a, b, x, l_ik, m)
    l = a - x(1:2) - R_z(x(3))*b;
    u = zeros(2,4);
    for i=1:m
        u(:,i) = l(:,i)/l_ik(i);
    end
end

function R = R_z(theta)
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
end