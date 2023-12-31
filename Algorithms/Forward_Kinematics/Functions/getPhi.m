function phi = getPhi(a, b, y, l, k_spring, m)
    phi = zeros(m,1);
    for i = 1:size(phi,1)
        phi(i) = nuSquared(a(:,i), b(:,i), y(:), l(i), k_spring);
    end
end

function nuSquared = nuSquared(a, b, y, l, k_spring)
    nuSquared = ...
        k_spring*( ...
              (a(1) - y(1) - b(1)*cos(y(3)) + b(2)*sin(y(3)))^2 ...
            + (y(2) - a(2) + b(2)*cos(y(3)) + b(1)*sin(y(3)))^2 ...
            - l^2 ...
        )^2;
end