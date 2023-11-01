% f = f_M + f_V
% f_M = 1/2*(f_min + f_max)


clc; clear all;

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


f_min = ones(4,1)*10; %[N]
f_max = ones(4,1)*25; %[N]
k_spring = 2000;
y = [-0.8; 1.1; 0.0];
l_ik = inverseKinematics(a,b,y,4);

u = cableUnitVectors(a,b,y,l_ik,4);
A_T = structureMatrix(u, y(3), b, 4);
A_T = A_T(2:4,:);
A_T_inv = pinv(A_T)
h = null(A_T);
w_p   = [0; 0; -9.81; 0.0; 0; 0];
w_p = w_p(2:4);

f_0 = A_T_inv*w_p;
lambda_l = max((f_min + f_0)./h);
lambda_h = min((f_max + f_0)./h);
sol_exists = lambda_l < lambda_h;
if sol_exists
    disp('Solution exists')
else
    disp('No solution exists')
end

f_M = 0.5*(f_min + f_max);
f_CFM = f_M - A_T_inv*(w_p + A_T*f_M)

i = -1;
for j=1:size(f_CFM)
    if f_CFM(j) < f_min(j)
        i = j;
        break;
    end
end

A_T
f_CFM
% A_T_marked is the structure matrix with the i-th column dropped
% f_CFM_marked is the force distribution vector with the i-th 
% element dropped
if i > 0
    A_T_marked = [A_T(:,1:i-1) A_T(:,i+1:end)];
%     f_CFM_marked = [f_CFM(1:i-1); f_CFM(i+1:end)];
    
    a_i = A_T(:,i);
    w_marked = f_min(i)*a_i + w_p;
    A_T_marked_inv = pinv(A_T_marked);
    f_M_marked = [f_M(1:i-1); f_M(i+1:end)];
    f_CFM_marked = f_M_marked - A_T_marked_inv*(w_marked + A_T_marked*f_M_marked)
end
%%
f_min = ones(4,1)*10; %[N]
f_max = ones(4,1)*25; %[N]
k_spring = 2000;
y = [-0.8; 1.1; 0.0];
l_ik = inverseKinematics(a,b,y,4);

u = cableUnitVectors(a,b,y,l_ik,4);
A_T = structureMatrix(u, y(3), b, 4);
A_T = A_T(2:4,:)
% A_T_inv = pinv(A_T)
% h = null(A_T);
w_p   = [0; 0; -9.81; 0.0; 0; 0];
w_p = w_p(2:4);
f_CFM = f_M - A_T_inv*(w_p + A_T*f_M)
f_CFM = improvedCFM(A_T, w_p, f_min, f_max, 0, 1)


function f_CFM = improvedCFM(A_T, w_p, f_min, f_max, r_level, max_r_level)
    r_level = r_level + 1;
    A_T_inv = pinv(A_T);
    f_M = 0.5*(f_min + f_max);
    f_CFM = f_M - A_T_inv*(w_p + A_T*f_M);
    % Is recursion level above maximum level?
    
    i = -1;
    j = -1;
    % Find first component that breaks minimum force requirement
    for k=1:size(f_CFM)
        if f_CFM(k) < f_min(k)
            i = k;
            break;
%         elseif f_CFM(k) > f_max(k)
%             j = k;
        end
    end
    
    if i > 0
        % No solution exists. Try to find new solution by moving force
        % component on nullspace of A^T

        % A_T_marked is the structure matrix with the i-th column dropped
        % f_CFM_marked is the force distribution vector with the i-th 
        % element dropped
        
        
        a_i = A_T(:,i);
        w_p_marked = f_min(i)*a_i + w_p;
%         disp('min,recursion level ')
%         disp(r_level)
%         f_min_marked = [f_min(1:i-1); f_min(i+1:end)]
%         f_max_marked = [f_max(1:i-1); f_max(i+1:end)]
        
        if i == 1
            A_T_marked = A_T(:,i+1:end);
            f_min_marked = f_min(i+1:end);
            f_max_marked = f_max(i+1:end);
        else
            A_T_marked = [A_T(:,1:i-1) A_T(:,i+1:end)];
            f_min_marked = [f_min(1:i-1); f_min(i+1:end)];
            f_max_marked = [f_max(1:i-1); f_max(i+1:end)];
        end
%         if r_level > max_r_level
%             disp('Max recursion level reached')
%             f_CFM = f_min(i);
%             return;
%         else
            f_CFM_marked = improvedCFM( ...
                                A_T_marked, ...
                                w_p_marked, ...
                                f_min_marked, ...
                                f_max_marked, ...
                                r_level, ...
                                max_r_level);
            f_CFM = [f_CFM_marked; f_min(i)];
%         end
    end
    for k=1:size(f_CFM)
        if f_CFM(k) > f_max(k)
            j = k;
            break;
        end
    end
    if j > 0
%         disp('max,recursion level ')
%         disp(r_level)
        if j == 1
            A_T_marked = A_T(:,j+1:end);
            f_min_marked = f_min(j+1:end);
            f_max_marked = f_max(j+1:end);
        else
            A_T_marked = [A_T(:,1:j-1) A_T(:,j+1:end)];
            f_min_marked = [f_min(1:j-1); f_min(j+1:end)];
            f_max_marked = [f_max(1:j-1); f_max(j+1:end)];
        end
        
        
        a_i = A_T(:,j);
        w_p_marked = f_min(j)*a_i + w_p;
%         if r_level > max_r_level
%             disp('Max recursion level reached')
%             f_CFM = f_max(i);
%             return;
%         else
            f_CFM_marked = improvedCFM( ...
                                A_T_marked, ...
                                w_p_marked, ...
                                f_min_marked, ...
                                f_max_marked, ...
                                r_level, ...
                                max_r_level);
            f_CFM = [f_CFM_marked; f_max(j)];
%         end
    end
    if r_level > max_r_level
        disp('Max recursion level reached')
%         f_CFM = f_min(i);
        return;
    end
end



% function f_CFM = improvedCFM(A_T, w_p, f_min, f_max, r_level, max_r_level)
%     r_level = r_level + 1;
%     A_T_inv = pinv(A_T);
%     f_M = 0.5*(f_min + f_max);
%     f_CFM = f_M - A_T_inv*(w_p + A_T*f_M);
%     % Is recursion level above maximum level?
%     
%     i = -1;
%     j = -1;
%     % Find first component that breaks minimum force requirement
%     for k=1:size(f_CFM)
%         if f_CFM(k) < f_min(k)
%             i = k;
%             break;
% %         elseif f_CFM(k) > f_max(k)
% %             j = k;
%         end
%     end
%     
%     if i > 0
%         % No solution exists. Try to find new solution by moving force
%         % component on nullspace of A^T
% 
%         % A_T_marked is the structure matrix with the i-th column dropped
%         % f_CFM_marked is the force distribution vector with the i-th 
%         % element dropped
%         
%         
%         a_i = A_T(:,i);
%         w_p_marked = f_min(i)*a_i + w_p;
% %         disp('min,recursion level ')
% %         disp(r_level)
% %         f_min_marked = [f_min(1:i-1); f_min(i+1:end)]
% %         f_max_marked = [f_max(1:i-1); f_max(i+1:end)]
%         
%         if i == 1
%             A_T_marked = A_T(:,i+1:end);
%             f_min_marked = f_min(i+1:end);
%             f_max_marked = f_max(i+1:end);
%         else
%             A_T_marked = [A_T(:,1:i-1) A_T(:,i+1:end)];
%             f_min_marked = [f_min(1:i-1); f_min(i+1:end)];
%             f_max_marked = [f_max(1:i-1); f_max(i+1:end)];
%         end
% %         if r_level > max_r_level
% %             disp('Max recursion level reached')
% %             f_CFM = f_min(i);
% %             return;
% %         else
%             f_CFM_marked = improvedCFM( ...
%                                 A_T_marked, ...
%                                 w_p_marked, ...
%                                 f_min_marked, ...
%                                 f_max_marked, ...
%                                 r_level, ...
%                                 max_r_level);
%             f_CFM = [f_CFM_marked; f_min(i)];
% %         end
%     end
%     for k=1:size(f_CFM)
%         if f_CFM(k) > f_max(k)
%             j = k;
%             break;
%         end
%     end
%     if j > 0
% %         disp('max,recursion level ')
% %         disp(r_level)
%         if j == 1
%             A_T_marked = A_T(:,j+1:end);
%             f_min_marked = f_min(j+1:end);
%             f_max_marked = f_max(j+1:end);
%         else
%             A_T_marked = [A_T(:,1:j-1) A_T(:,j+1:end)];
%             f_min_marked = [f_min(1:j-1); f_min(j+1:end)];
%             f_max_marked = [f_max(1:j-1); f_max(j+1:end)];
%         end
%         
%         
%         a_i = A_T(:,j);
%         w_p_marked = f_min(j)*a_i + w_p;
% %         if r_level > max_r_level
% %             disp('Max recursion level reached')
% %             f_CFM = f_max(i);
% %             return;
% %         else
%             f_CFM_marked = improvedCFM( ...
%                                 A_T_marked, ...
%                                 w_p_marked, ...
%                                 f_min_marked, ...
%                                 f_max_marked, ...
%                                 r_level, ...
%                                 max_r_level);
%             f_CFM = [f_CFM_marked; f_max(j)];
% %         end
%     end
%     if r_level > max_r_level
%         disp('Max recursion level reached')
% %         f_CFM = f_min(i);
%         return;
%     end
% end




% function f_CFM = improvedCFM(A_T, w_p, f_min, f_max)
%     A_T_inv = pinv(A_T)
%     h = null(A_T)
%     f_M = 0.5*(f_min + f_max);
%     
%     f_CFM = f_M - A_T_inv*(w_p + A_T*f_M);
% 
%     % Check for feasible solutions:
%     f_0 = A_T_inv*w_p;
%     lambda_l = max((f_min + f_0)./h);
%     lambda_h = min((f_max + f_0)./h);
%     if lambda_l < lambda_h
%         disp('Solution exists')
%         % Solution exists. Return force distribution vector
%         
%         return;
%     else
%         disp('No solution exists')
%         % No solution exists. Try to find new solution by moving force
%         % component on nullspace of A^T
%         i = -1;
%         % Find first component that breaks minimum force requirement
%         for j=1:size(f_CFM)
%             if f_CFM(j) < f_min(j)
%                 i = j;
%                 break;
%             end
%         end
%         % A_T_marked is the structure matrix with the i-th column dropped
%         % f_CFM_marked is the force distribution vector with the i-th 
%         % element dropped
%         if i > 0
%             A_T_marked = [A_T(:,1:i-1) A_T(:,i+1:end)];
%         %     f_CFM_marked = [f_CFM(1:i-1); f_CFM(i+1:end)];
%             
%             a_i = A_T(:,i);
%             w_p_marked = f_min(i)*a_i + w_p;
% %             A_T_marked_inv = pinv(A_T_marked);
% %             f_M_marked = [f_M(1:i-1); f_M(i+1:end)];
%             
% %             f_CFM_marked = f_M_marked - A_T_marked_inv*(w_p_marked + A_T_marked*f_M_marked)
%             f_CFM = improvedCFM(A_T_marked, w_p_marked, f_min, f_max)
%             disp('hey')
% %             f_M_marked - A_T_marked_inv*(w_p_marked + A_T_marked*f_M_marked)
%         end
%     end
% end