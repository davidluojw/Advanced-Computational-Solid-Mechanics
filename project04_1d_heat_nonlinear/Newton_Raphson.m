% This is the Newton-Raphson algorithm solver for HW3
clear all; clc;

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


for n = 0:load_num
    iter_num = 0;    % record the iterate numbers

    R(1) = R1(d(1),d(2),n);
    R(2) = R2(d(1),d(2),n);
    R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
    R_norm0 = R_norm;

    while true
        J(1,1) = dN1_dd1(d(1),d(2));
        J(1,2) = dN1_dd2(d(1),d(2));
        J(2,1) = dN2_dd1(d(1),d(2));
        J(2,2) = dN2_dd2(d(1),d(2));

        Deltad = -J\R;

        d = d + Deltad;
        iter_num = iter_num + 1;

        if iter_num == 15
            break;
        end

        R(1) = R1(d(1),d(2),n);
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));

        if R_norm <= tol * R_norm0
            break;
        end
    end
    d_sol(:,n+1) = d;
end


% plot (a) x =15
figure;
d_exact = linspace(-1,7,100);
N1_exact = N1(d_exact, d_exact);
plot(d_exact, N1_exact, LineWidth=2);
title("N1 vs d1");
xlabel("d1");
ylabel("N1");
legend("exact, x=15",'Location', 'Best', 'FontSize', 14, 'Box', 'on');

%plot (b) x=15
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

% plot (a) x = 25
figure;
x = 25;
N1 = @(d1,d2) x * d1 ./ (10.0 - d1) - 0.5 * d2.*d2;
N2 = @(d1,d2) d2 - d1;
N1_exact = N1(d_exact, d_exact);
plot(d_exact, N1_exact, LineWidth=2);
title("N1 vs d1");
xlabel("d1");
ylabel("N1");
legend("exact, x=25",'Location', 'Best', 'FontSize', 14, 'Box', 'on');


%plot (b) x=25
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







