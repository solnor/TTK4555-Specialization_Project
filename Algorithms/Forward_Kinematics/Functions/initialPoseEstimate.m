function y_0 = initialPoseEstimate(a, b, l)
    m = size(a,2);
    r_low_i  = zeros(2, m);
    r_high_i = zeros(2, m);
    
    for i=1:m
        r_low_i(:, i)  = ...
            a(:,i) - (l(i) + norm(b(:,i), 2))*[1;1];
        r_high_i(:, i) = ...
            a(:,i) + (l(i) + norm(b(:,i), 2))*[1;1];
    end
    
    r_low = [max(r_low_i(1,:)); 
             max(r_low_i(2,:))];
    
    r_high = [min(r_high_i(1,:));
              min(r_high_i(2,:));];
    
    r_0 = [0.5*(r_low(1) + r_high(1)); 
           0.5*(r_low(2) + r_high(2))];
    y_0 = [r_0; 0.0];
end