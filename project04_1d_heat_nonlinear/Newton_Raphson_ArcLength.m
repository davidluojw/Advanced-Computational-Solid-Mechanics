
% This is the Newton-Raphson algorithm solver for HW3
clear all; clc;

%% x =10
% exact internal force
x = 10;
N1 = @(d1,d2) x * d1 ./ (10.0 - d1) - 0.5 * d2.*d2;
N2 = @(d1,d2) d2 - d1;

% Jacobian
J = zeros(2,2);
dN1_dd1 = @(d1,d2) 10*x / ((10.0-d1)*(10.0-d1));
dN1_dd2 = @(d1,d2) -d2;
dN2_dd1 = @(d1,d2) -1.0;
dN2_dd2 = @(d1,d2) 1.0;


% Arc length
point_num = 40;     % point numbers
delta_a = 0.3;      % arc-length increment

% arc-length function
% deltad: the total change over the step to current iterate

arclength_fun = @(deltad, deltalambda) deltad' *  deltad + deltalambda * deltalambda - delta_a * delta_a;


% external force with arc-length method
load_max = 10.0;
F1 = @() load_max;
F2 = @() 0.0;
F_ext =[F1(); F2()];

% Residual R with arc-length method
R = zeros(2,1);
R1 = @(d1,d2,lambda) lambda * F1() - N1(d1,d2);
R2 = @(d1,d2,lambda) lambda * F2() - N2(d1,d2);


% Newton-Raphson Algorithm with using arc-length method

% increment of N-R iteration
Deltad = zeros(2,1);
Deltalambda = 0;

% increment of the load step
deltad = zeros(2,1);
deltalambda = 0;

% initial point value
d = zeros(2,1);
lambda = 0.0;

% tolerance
tol = 1.0e-8;

%  solutions recording 
d_sol = zeros(2, point_num + 1);            % stored convergent d_n+1 at every load step
d_sol(:,1) = d;
lambda_sol = zeros(1, point_num + 1);      % stored convergent lambda_n+1 at every load step
lambda_sol(:, 1) = lambda;

deltad_set = zeros(2, point_num + 1);       % stored deltad_n^(0)) at every load step
deltad_set(:,1) = deltad;
deltalambda_set = zeros(1, point_num + 1);  % stored deltalambda_n^(0) at every load step
deltalambda_set(:, 1) = deltalambda;  
iter_set = zeros(1, point_num+1);     % stored convergent iterate number at every load step
detJ_set = zeros(1, point_num + 1);   % stored detK(d_n^(0)) at every load step


% Arc-length process
for n = 1:point_num

    % Initialization, Predictor
    d_n = d;
    lambda_n = lambda;

    % K(d_n+1^(0))
    J(1,1) = dN1_dd1(d(1),d(2));
    J(1,2) = dN1_dd2(d(1),d(2));
    J(2,1) = dN2_dd1(d(1),d(2));
    J(2,2) = dN2_dd2(d(1),d(2));
    K0 = J;
    detJ = det(J);   % det K(d_n+1^(0))
    detJ_set(:, n) = detJ;
    
     
    q_bar = K0 \ F_ext;          % q_bar = K(d_n+1^(0))^(-1) * F_ext    

    if n == 1
        % 手动设置
        Deltalambda = sqrt( delta_a * delta_a / (1 + q_bar' * q_bar));  % Deltalambda^(0)
    else
        % the sign of Deltalambda is based on the sign of detK at previous
        % load step 
        if sign(detJ) == -sign(detJ_set(:,n-1))
            Deltalambda = -sign(deltalambda_set(:, n-1)) ... 
                * sqrt( delta_a * delta_a / (1 + q_bar' * q_bar));   % Deltalambda^(0)
        else
            Deltalambda = sign(deltalambda_set(:, n-1)) ... 
                * sqrt( delta_a * delta_a / (1 + q_bar' * q_bar));    % Deltalambda^(0)
        end
    end

    Deltad = Deltalambda * q_bar;  % Deltad^(0) = Deltalambda^(0)* q_bar

    % Update points
    d = d + Deltad;   % d_n+1^(1) = d_n+1^(0) + Deltad^(0)
    lambda = lambda + Deltalambda;   % lambda_n+1^(1) = lambda_n+1^(0) + Deltalambda^(0)

    deltad = Deltad;              % deltad^(0)
    deltalambda = Deltalambda;    % deltalambda^(0)
    deltad_set(:, n) = deltad;
    deltalambda_set(:,n) = deltalambda; 
 

    iter_step = 1;   % iteration step
    % N-R iteration with constraint surface
    while true

        q_bar = K0 \ F_ext;     % q_bar = K(d_n+1^(0))^(-1) * F_ext 

        % R(d_n+1^(i), lambda_n+1^(i))
        R(1) = R1(d(1),d(2),lambda);
        R(2) = R2(d(1),d(2),lambda);

        Deltad_bar = K0 \ R;    % Deltad_bar = K(d_n+1^(0))^(-1) * R(d_n+1^(i), lambda_n+1^(i))

        deltad = d - d_n;                 % deltad^(i) = d_n+1^(i) - d_n
        deltalambda = lambda - lambda_n;  % deltalambda^(i) = lambda_n+1^(i) - lambda_n

        f_np1 = arclength_fun(deltad, deltalambda);  % f(d_n+1^(i), lambda_n+1^(0))

        % Deltalambda^(i)
        Deltalambda = ( - f_np1 - 2 * deltad' * Deltad_bar ) / ( 2 * deltad' * q_bar + 2 * deltalambda * deltalambda);

        % Deltad^(i) = Deltalambda^(i) * q_bar + Deltad_bar
        Deltad = Deltalambda * q_bar + Deltad_bar;

        % Update points
        d = d + Deltad;                 % d_n+1^(i+1) = d_n+1^(i) + Deltad^(i)
        lambda = lambda + Deltalambda;  % lambda_n+1^(i+1) = lambda_n+1^(i) + Deltalambda^(i)

        iter_step = iter_step + 1;
    
        if iter_step == 100
            break;
        end

        % R(d_n+1^(i+1), lambda_n+1^(i+1)) = R_n+1^(i+1)
        R(1) = R1(d(1),d(2),lambda);
        R(2) = R2(d(1),d(2),lambda);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
    
        if R_norm <= tol
            break;
        end

    end

    iter_set(:,n) = iter_step;

    % *: convergent
    d_sol(:,n) = d;                      % d_n+1*
    lambda_sol(:,n) = lambda;            % lambda_n+1*
end



% plot (a) x = 10
figure;
d_exact = linspace(-1,8.5,100);
N1_exact = N1(d_exact, d_exact);
plot(d_exact, N1_exact, LineWidth=2);
title("N1 vs d1");
xlabel("d1");
ylabel("N1");
legend("exact, x=10",'Location', 'Best', 'FontSize', 14, 'Box', 'on');

%plot (b) x = 10
figure;
d_exact = linspace(-1,8.5,100);
N1_exact = N1(d_exact, d_exact);
N1_sol = N1(d_sol(1,:), d_sol(1,:));
plot(d_exact, N1_exact, LineWidth=2);
hold on;
scatter(d_sol(1,:),N1_sol,"filled");
title("N1 vs d1 (x=10)");
xlabel("d1");
ylabel("N1");
legend("exact","numerical",'Location', 'Best', 'FontSize', 14, 'Box', 'on');

%plot (c)  x= 10
figure;
load_steps = linspace(1,point_num + 1, point_num + 1);
scatter(load_steps, iter_set, "filled");
title("iterations vs load step (x=10)");
% xlim([0,]);
% ylim([0,16]);
xlabel("load step");
ylabel("iterations");
legend("x = 10",'Location', 'Best', 'FontSize', 14, 'Box', 'on');


