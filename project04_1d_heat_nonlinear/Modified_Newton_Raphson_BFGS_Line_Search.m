% This is the Newton-Raphson algorithm solver for HW3
clear all; clc;

%% x =15
% exact internal force
x = 15;
N1 = @(d1,d2) x * d1 ./ (10.0 - d1) - 0.5 * d2.*d2;
N2 = @(d1,d2) d2 - d1;

% external force
load_step = 0.25;
load_max = 10.0;
load_num = load_max / load_step;
F1 = @(n) load_step * n;
F2 = @(n) 0.0;

% Jacobian
J = zeros(2,2);
dN1_dd1 = @(d1,d2) 10*x / ((10.0-d1)*(10.0-d1));
dN1_dd2 = @(d1,d2) -d2;
dN2_dd1 = @(d1,d2) -1.0;
dN2_dd2 = @(d1,d2) 1.0;

% Residual R
R = zeros(2,1);
R1 = @(d1,d2,n) N1(d1,d2) - F1(n);
R2 = @(d1,d2,n) N2(d1,d2) - F2(n);


% Newton-Raphson Algorithm
% K * Deltad = -R
Deltad = zeros(2,1);
d = zeros(2,1);
tol = 1.0e-4;
d_sol = zeros(2, load_num + 1);
d_sol(:,1) = d;
iter_num = zeros(1, load_num+1);
max_iter = 15;

% update vectors: V and W
V = zeros(2, max_iter + 1);
W = zeros(2, max_iter + 1);

% line search
stol = 0.5;


for n = 1:load_num
    iter_step = 0;    % record the iterate numbers

    % get the difference between first two residual
    % Delta R = R(i+1) - R(i)
    R(1) = R1(d(1),d(2),n);
    R(2) = R2(d(1),d(2),n);
    R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
    R_norm0 = R_norm;

    R_old = R;   % R(i)
    d_old = d;   % d(i)

    J(1,1) = dN1_dd1(d(1),d(2));
    J(1,2) = dN1_dd2(d(1),d(2));
    J(2,1) = dN2_dd1(d(1),d(2));
    J(2,2) = dN2_dd2(d(1),d(2));

    Deltad = -J\R;   % Deltad(i)
    d = d + Deltad;  % d(i+1) = d(i) + Deltad(i)

    % this is the new residual R(i+1)
    R(1) = R1(d(1),d(2),n);
    R(2) = R2(d(1),d(2),n);

    % BFGS start, solve K_bar * Delta d(i+1) = R(i+1)
    % where Delta d(i+1) = d(i+2) - d(i+1)
    % actually, BFGS is to find d(i+2)
    while true
        iter_step = iter_step + 1;

        G0 = Deltad' * R;

        R_tilde = R;  % R_tilde(0) = R(i+1)

        % compute update vectors V and W
        % V(i-k+1) = Deltad(i-k+1) / (Deltad(i-k+1) * DeltaR(i-k+1))
        % DeltaR(i-k+1) = R(i-k+2) - R(i-k+1)
        % Deltad(i-k+1) = d(i-k+2) - d(i-k+1)

        % k = 1, i-k+1 = i, i-k+2 = i+1
        DeltaR_k = R_tilde - R_old;     % Delta R(i-k+1) = Delta R(i) = R(i+1) - R(i)
        Deltad_k = Deltad;              % Delta d(i-k+1) = Delta d(i) 

        % V(i-k+1) = Deltad(i-k+1) / (Deltad(i-k+1) * DeltaR(i-k+1))
        V(:,iter_step) = Deltad_k / (Deltad_k' * DeltaR_k);

        % alpha = (-Delta R(i-k+1) * Delta d(i-k+1) / (R(i-k+1) * Delta d(i-k+1)) )^(1/2)
        %       = (-Delta R(i) * Delta d(i) / (R(i) * Delta d(i) ) )^(1/2)
        alpha = sqrt(-DeltaR_k' * Deltad_k / (R_old' * Deltad_k));

        % W(i-k+1) = alpha * R(i-k+1) - Delta R(i-k+1)
        %          = alpha * R(i) - Delta R(i)
        W(:,iter_step) = alpha * R_old - DeltaR_k;

        % Right-side updates
        for k = 1:iter_step
            R_tilde = R_tilde + (V(:, iter_step - k + 1)'* R_tilde) * W(:,iter_step - k + 1);
        end

        % Solve intermediate system
        % Deltad_tilde(0) = J^(-1) * R_tilde(i)
        Deltad_tilde = -J\R_tilde;

        % Left-side update
        for k = 1:iter_step
            Deltad_tilde = Deltad_tilde + (W(:, k)' * Deltad_tilde) * V(:, k);
        end
        d_old = d;
        d = d + Deltad_tilde;  % d(i+2) = d(i+1) + Deltad_tilde
        Deltad = d - d_old;    % Deltad(i+1) = d(i+2) - d(i+1)

        % R(i+2)
        R(1) = R1(d(1),d(2),n);  
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));

        G = Deltad' * R;

        % Line Search
        s = LineSearch(G0,G,d,R1,R2,Deltad,n,stol);

        % update the new updates
        d = d - Deltad;
        d = d + s * Deltad;
    
        R(1) = R1(d(1),d(2),n);
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));

        if R_norm <= tol * R_norm0
            break;
        end
        if iter_step == max_iter
            break;
        end
    end  % end of iteration
    iter_num(:,n+1) = iter_step;
    d_sol(:,n+1) = d;
end  % end of load step


% plot (a) x = 15
figure;
d_exact = linspace(-1,7,100);
N1_exact = N1(d_exact, d_exact);
plot(d_exact, N1_exact, LineWidth=2);
title("N1 vs d1");
xlabel("d1");
ylabel("N1");
legend("exact, x=15",'Location', 'Best', 'FontSize', 14, 'Box', 'on');

%plot (b) x = 15
figure;
d_exact = linspace(-1,7,100);
N1_exact = N1(d_exact, d_exact);
N1_sol = N1(d_sol(1,:), d_sol(1,:));
plot(d_exact, N1_exact, LineWidth=2);
hold on;
scatter(d_sol(1,:),N1_sol,"filled");
title("N1 vs d1 (x=15)");
xlabel("d1");
ylabel("N1");
legend("exact","numerical",'Location', 'Best', 'FontSize', 14, 'Box', 'on');

%plot (c)  x= 15
figure;
load_steps = linspace(0,load_num + 1, load_num + 1);
scatter(load_steps, iter_num, "filled");
title("iterations vs load_step (x=15)");
xlim([0,42]);
ylim([0,16]);
xlabel("load_step");
ylabel("iterations");
legend("x = 15",'Location', 'Best', 'FontSize', 14, 'Box', 'on');


%% x = 25
% exact internal force
x = 25;
N1 = @(d1,d2) x * d1 ./ (10.0 - d1) - 0.5 * d2.*d2;
N2 = @(d1,d2) d2 - d1;

% external force
load_step = 0.25;
load_max = 10.0;
load_num = load_max / load_step;
F1 = @(n) load_step * n;
F2 = @(n) 0.0;

% Jacobian
J = zeros(2,2);
dN1_dd1 = @(d1,d2) 10*x / ((10.0-d1)*(10.0-d1));
dN1_dd2 = @(d1,d2) -d2;
dN2_dd1 = @(d1,d2) -1.0;
dN2_dd2 = @(d1,d2) 1.0;

% Residual R
R = zeros(2,1);
R1 = @(d1,d2,n) N1(d1,d2) - F1(n);
R2 = @(d1,d2,n) N2(d1,d2) - F2(n);


% Newton-Raphson Algorithm
% K * Deltad = -R
Deltad = zeros(2,1);
d = zeros(2,1);
tol = 1.0e-4;
d_sol = zeros(2, load_num + 1);
d_sol(:,1) = d;
iter_num = zeros(1, load_num+1);
max_iter = 15;

% update vectors: V and W
V = zeros(2, max_iter + 1);
W = zeros(2, max_iter + 1);

% line search
stol = 0.5;


for n = 1:load_num
    iter_step = 0;    % record the iterate numbers

    % get the difference between first two residual
    % Delta R = R(i+1) - R(i)
    R(1) = R1(d(1),d(2),n);
    R(2) = R2(d(1),d(2),n);
    R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
    R_norm0 = R_norm;

    R_old = R;   % R(i)
    d_old = d;   % d(i)

    J(1,1) = dN1_dd1(d(1),d(2));
    J(1,2) = dN1_dd2(d(1),d(2));
    J(2,1) = dN2_dd1(d(1),d(2));
    J(2,2) = dN2_dd2(d(1),d(2));

    Deltad = -J\R;   % Deltad(i)
    d = d + Deltad;  % d(i+1) = d(i) + Deltad(i)

    % this is the new residual R(i+1)
    R(1) = R1(d(1),d(2),n);
    R(2) = R2(d(1),d(2),n);

    % BFGS start, solve K_bar * Delta d(i+1) = R(i+1)
    % where Delta d(i+1) = d(i+2) - d(i+1)
    % actually, BFGS is to find d(i+2)
    while true
        iter_step = iter_step + 1;

        G0 = Deltad' * R;

        R_tilde = R;  % R_tilde(0) = R(i+1)

        % compute update vectors V and W
        % V(i-k+1) = Deltad(i-k+1) / (Deltad(i-k+1) * DeltaR(i-k+1))
        % DeltaR(i-k+1) = R(i-k+2) - R(i-k+1)
        % Deltad(i-k+1) = d(i-k+2) - d(i-k+1)

        % k = 1, i-k+1 = i, i-k+2 = i+1
        DeltaR_k = R_tilde - R_old;     % Delta R(i-k+1) = Delta R(i) = R(i+1) - R(i)
        Deltad_k = Deltad;              % Delta d(i-k+1) = Delta d(i) 

        % V(i-k+1) = Deltad(i-k+1) / (Deltad(i-k+1) * DeltaR(i-k+1))
        V(:,iter_step) = Deltad_k / (Deltad_k' * DeltaR_k);

        % alpha = (-Delta R(i-k+1) * Delta d(i-k+1) / (R(i-k+1) * Delta d(i-k+1)) )^(1/2)
        %       = (-Delta R(i) * Delta d(i) / (R(i) * Delta d(i) ) )^(1/2)
        alpha = sqrt(-DeltaR_k' * Deltad_k / (R_old' * Deltad_k));

        % W(i-k+1) = alpha * R(i-k+1) - Delta R(i-k+1)
        %          = alpha * R(i) - Delta R(i)
        W(:,iter_step) = alpha * R_old - DeltaR_k;

        % Right-side updates
        for k = 1:iter_step
            R_tilde = R_tilde + (V(:, iter_step - k + 1)'* R_tilde) * W(:,iter_step - k + 1);
        end

        % Solve intermediate system
        % Deltad_tilde(0) = J^(-1) * R_tilde(i)
        Deltad_tilde = -J\R_tilde;

        % Left-side update
        for k = 1:iter_step
            Deltad_tilde = Deltad_tilde + (W(:, k)' * Deltad_tilde) * V(:, k);
        end
        d_old = d;
        d = d + Deltad_tilde;  % d(i+2) = d(i+1) + Deltad_tilde
        Deltad = d - d_old;    % Deltad(i+1) = d(i+2) - d(i+1)

        % R(i+2)
        R(1) = R1(d(1),d(2),n);  
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));

        G = Deltad' * R;

        % Line Search
        s = LineSearch(G0,G,d,R1,R2,Deltad,n,stol);

        % update the new updates
        d = d - Deltad;
        d = d + s * Deltad;
    
        R(1) = R1(d(1),d(2),n);
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));

        if R_norm <= tol * R_norm0
            break;
        end
        if iter_step == max_iter
            break;
        end
    end  % end of iteration
    iter_num(:,n+1) = iter_step;
    d_sol(:,n+1) = d;
end  % end of load step


% plot (a) x = 25
figure;
d_exact = linspace(-1,7,100);
N1_exact = N1(d_exact, d_exact);
plot(d_exact, N1_exact, LineWidth=2);
title("N1 vs d1");
xlabel("d1");
ylabel("N1");
legend("exact, x=25",'Location', 'Best', 'FontSize', 14, 'Box', 'on');

%plot (b) x = 25
figure;
d_exact = linspace(-1,7,100);
N1_exact = N1(d_exact, d_exact);
N1_sol = N1(d_sol(1,:), d_sol(1,:));
plot(d_exact, N1_exact, LineWidth=2);
hold on;
scatter(d_sol(1,:),N1_sol,"filled");
title("N1 vs d1 (x=25)");
xlabel("d1");
ylabel("N1");
legend("exact","numerical",'Location', 'Best', 'FontSize', 14, 'Box', 'on');

%plot (c)  x= 25
figure;
load_steps = linspace(0,load_num + 1, load_num + 1);
scatter(load_steps, iter_num, "filled");
title("iterations vs load_step (x=15)");
xlim([0,42]);
ylim([0,16]);
xlabel("load_step");
ylabel("iterations");
legend("x = 25",'Location', 'Best', 'FontSize', 14, 'Box', 'on');