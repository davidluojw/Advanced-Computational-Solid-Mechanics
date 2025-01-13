% This is the Newton-Raphson algorithm solver for HW3
clear all; clc;

% exact internal force
x = 15;
N1 = @(d1,d2) x * d1 / (10.0 - d1) - 0.5 * d2*d2;
N2 = @(d1,d2) d2 - d1;
load_step = 0.25;
load_max = 10.0;
load_num = load_max / load_step;


% external force
F1 = @(n) 0.25 * n;
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
iter_num = 0;

for n = 0:load_num
    R(1) = R1(d(1),d(2),n);
    R(2) = R2(d(1),d(2),n);
    R_norm = sqrt(R1 * R1 + R2 * R2);
    R_norm0 = R_norm;

    while true
        J(1,1) = dN1dd1(d(1),d(2));
        J(1,2) = dN1dd2(d(1),d(2));
        J(2,1) = dN2dd1(d(1),d(2));
        J(2,2) = dN2dd2(d(1),d(2));

        Deltad = -J\r;

        d = d + Deltad;
        iter_num = iter_num + 1;

        if iter_num == 15
            break;
        end

        R(1) = R1(d(1),d(2),n);
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R1 * R1 + R2 * R2);

        if R_norm <= tol * R_norm0
            break;
        end

    end
end






