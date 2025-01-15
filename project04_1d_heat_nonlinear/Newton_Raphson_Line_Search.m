% This is the Newton-Raphson with line search algorithm solver for HW3
% Reference: H. Matthies and G. Strange. IJNME, 14 (1979) 1613-1626
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

% line search
stol = 0.5;


for n = 0:load_num
    iter_step = 0;    % record the iterate numbers

    R(1) = R1(d(1),d(2),n);
    R(2) = R2(d(1),d(2),n);
    R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
    R_norm0 = R_norm;

    while true
        J(1,1) = dN1_dd1(d(1),d(2));
        J(1,2) = dN1_dd2(d(1),d(2));
        J(2,1) = dN2_dd1(d(1),d(2));
        J(2,2) = dN2_dd2(d(1),d(2));

        Deltad = -J\R;     % Deltad(i)
        d_temp = d;        % d_temp = d(i)
        G0 = Deltad' * R;   % G0 = Deltad(i) * R(i)

        d = d + Deltad;     % d(i+1) = d(i) + Deltad(i)
        iter_step = iter_step + 1;

        if iter_step == 15
            break;
        end

        % R(i+1)
        R(1) = R1(d(1),d(2),n);
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));

        G = Deltad' * R;    % Deltad(i) * R(i+1) 

        % Line Search
        s = LineSearch(G0,G,d_temp,R1,R2,Deltad,n,stol);  % d(i), Deltad(i)

        % update the new updates
        d = d_temp + s * Deltad;   % d(i+1) = d(i) + s*Deltad(i)

        Deltad = d - d_temp;       % Delta(i) = d(i+1) - d(i)

        % R(i+1)
        R(1) = R1(d(1),d(2),n);
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));

        if R_norm <= tol * R_norm0
            break;
        end
    end
    iter_num(:,n+1) = iter_step;
    d_sol(:,n+1) = d;
end


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

% line search
stol = 0.5;


for n = 0:load_num
    iter_step = 0;    % record the iterate numbers

    % R(i)
    R(1) = R1(d(1),d(2),n);
    R(2) = R2(d(1),d(2),n);
    R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
    R_norm0 = R_norm;

    while true
        J(1,1) = dN1_dd1(d(1),d(2));
        J(1,2) = dN1_dd2(d(1),d(2));
        J(2,1) = dN2_dd1(d(1),d(2));
        J(2,2) = dN2_dd2(d(1),d(2));

        Deltad = -J\R;     % Deltad(i)
        d_temp = d;        % d_temp = d(i)
        G0 = Deltad' * R;   % G0 = Deltad(i) * R(i)

        d = d + Deltad;     % d(i+1) = d(i) + Deltad(i)
        iter_step = iter_step + 1;

        if iter_step == 15
            break;
        end

        % R(i+1)
        R(1) = R1(d(1),d(2),n);
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));

        G = Deltad' * R;    % Deltad(i) * R(i+1) 

        % Line Search
        s = LineSearch(G0,G,d_temp,R1,R2,Deltad,n,stol);  % d(i), Deltad(i)

        % update the new updates
        d = d_temp + s * Deltad;   % d(i+1) = d(i) + s*Deltad(i)

        Deltad = d - d_temp;       % Delta(i) = d(i+1) - d(i)

        % R(i+1)
        R(1) = R1(d(1),d(2),n);
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));


        if R_norm <= tol * R_norm0
            break;
        end
    end
    iter_num(:,n+1) = iter_step;
    d_sol(:,n+1) = d;
end


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
title("iterations vs load_step (x=25)");
xlim([0,42]);
ylim([0,16]);
xlabel("load_step");
ylabel("iterations");
legend("x = 25",'Location', 'Best', 'FontSize', 14, 'Box', 'on');