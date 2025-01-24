% This is the Newton-Raphson algorithm solver for HW3
clear all; clc;
% 
% %% x =10
% % exact internal force
% x = 10;
% N1 = @(d1,d2) x * d1 ./ (10.0 - d1) - 0.5 * d2.*d2;
% N2 = @(d1,d2) d2 - d1;
% 
% % external force
% load_step = 0.25;
% load_max = 10.0;
% load_num = load_max / load_step;
% F1 = @(n) load_step * n;
% F2 = @(n) 0.0;
% 
% % Jacobian
% J = zeros(2,2);
% dN1_dd1 = @(d1,d2) 10*x / ((10.0-d1)*(10.0-d1));
% dN1_dd2 = @(d1,d2) -d2;
% dN2_dd1 = @(d1,d2) -1.0;
% dN2_dd2 = @(d1,d2) 1.0;
% 
% % Residual R
% R = zeros(2,1);
% R1 = @(d1,d2,n) N1(d1,d2) - F1(n);
% R2 = @(d1,d2,n) N2(d1,d2) - F2(n);
% 
% 
% % Newton-Raphson Algorithm with using arc-length method
% % [ K , -F^ext; (df/d delta_d)^T, df/d delta_lambda ] * 
% % [ Deltad; Deltalambda] = [lambda * F^ext - F^int(d); delta_a - f ]
% % f(delta_d, delta_lambda) = ()
% Deltad = zeros(2,1);
% d = zeros(2,1);
% tol = 1.0e-4;
% d_sol = zeros(2, load_num + 1);
% d_sol(:,1) = d;
% iter_num = zeros(1, load_num+1);
% 
% 
% for n = 0:load_num
%     iter_step = 0;    % record the iterate numbers
% 
%     R(1) = R1(d(1),d(2),n);
%     R(2) = R2(d(1),d(2),n);
%     R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
%     R_norm0 = R_norm;
% 
%     while true
%         J(1,1) = dN1_dd1(d(1),d(2));
%         J(1,2) = dN1_dd2(d(1),d(2));
%         J(2,1) = dN2_dd1(d(1),d(2));
%         J(2,2) = dN2_dd2(d(1),d(2));
% 
%         Deltad = -J\R;
% 
%         d = d + Deltad;
%         iter_step = iter_step + 1;
% 
%         if iter_step == 15
%             break;
%         end
% 
%         R(1) = R1(d(1),d(2),n);
%         R(2) = R2(d(1),d(2),n);
%         R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
% 
%         if R_norm <= tol * R_norm0
%             break;
%         end
%     end
%     iter_num(:,n+1) = iter_step;
%     d_sol(:,n+1) = d;
% end
% 
% 
% % plot (a) x = 10
% figure;
% d_exact = linspace(-1,8.5,100);
% N1_exact = N1(d_exact, d_exact);
% plot(d_exact, N1_exact, LineWidth=2);
% title("N1 vs d1");
% xlabel("d1");
% ylabel("N1");
% legend("exact, x=10",'Location', 'Best', 'FontSize', 14, 'Box', 'on');
% 
% %plot (b) x = 10
% figure;
% d_exact = linspace(-1,8.5,100);
% N1_exact = N1(d_exact, d_exact);
% N1_sol = N1(d_sol(1,:), d_sol(1,:));
% plot(d_exact, N1_exact, LineWidth=2);
% hold on;
% scatter(d_sol(1,:),N1_sol,"filled");
% title("N1 vs d1 (x=10)");
% xlabel("d1");
% ylabel("N1");
% legend("exact","numerical",'Location', 'Best', 'FontSize', 14, 'Box', 'on');
% 
% %plot (c)  x= 10
% figure;
% load_steps = linspace(0,load_num + 1, load_num + 1);
% scatter(load_steps, iter_num, "filled");
% title("iterations vs load_step (x=10)");
% xlim([0,42]);
% ylim([0,16]);
% xlabel("load_step");
% ylabel("iterations");
% legend("x = 10",'Location', 'Best', 'FontSize', 14, 'Box', 'on');



% 定义初始条件
d = [0; 0]; % 初始位移
F_ext = [0; 0]; % 初始外部力
x = 1; % 参数x的值
arc_length = 0.1; % 弧长步长

% 定义增量加载步长
load_steps = 0:0.25:10;

% 迭代求解
for F1_ext = load_steps
    F_ext(1) = F1_ext;
    d = solve_nonlinear_system_with_arc_length(d, F_ext, x, arc_length);
    disp(['Load step: ', num2str(F1_ext), ', Solution: ', num2str(d')]);
end

function d = solve_nonlinear_system_with_arc_length(d, F_ext, x, arc_length)
    % 使用弧长法和牛顿迭代法求解非线性方程组
    tol = 1e-6;
    max_iter = 100;
    lambda = 0; % 初始化弧长参数

    for iter = 1:max_iter
        N = nonlinear_equations(d, x);
        J = jacobian_matrix(d, x);
        
        % 计算弧长约束
        arc_constraint = norm([d; lambda]) - arc_length;
        
        % 扩展雅可比矩阵以包含弧长约束
        J_ext = [J, zeros(2, 1); 2*d', 2*lambda];
        
        % 扩展非线性方程组以包含弧长约束
        N_ext = [F_ext - N; arc_constraint];
        
        % 计算增量
        delta = J_ext \ N_ext;
        d = d + delta(1:2);
        lambda = lambda + delta(3);
        
        if norm(delta) < tol
            break;
        end
    end
end

function N = nonlinear_equations(d, x)
    d1 = d(1);
    d2 = d(2);
    N1 = (x * d1) / (10 - d1) - 0.5 * d2^2;
    N2 = d2 - d1;
    N = [N1; N2];
end

function J = jacobian_matrix(d, x)
    d1 = d(1);
    d2 = d(2);
    J11 = (10 * x) / (10 - d1)^2;
    J12 = -d2;
    J21 = -1;
    J22 = 1;
    J = [J11, J12; J21, J22];
end