function l = inverseKinematics(a, b, pose, m)
    x     = pose(1);
    y     = pose(2);
    theta = pose(3);

    l = zeros(m,1);
    for i=1:m
        l(i) = norm(a(:,i)-[x;y]-R_z(theta)*b(:,i),2);
    end
end


function R = R_z(theta)
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
end