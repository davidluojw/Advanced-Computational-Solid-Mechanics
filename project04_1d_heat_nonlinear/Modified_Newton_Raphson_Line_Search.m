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

% line search
stol = 0.5;


for n = 0:load_num
    iter_step = 0;    % record the iterate numbers

    R(1) = R1(d(1),d(2),n);
    R(2) = R2(d(1),d(2),n);
    R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
    R_norm0 = R_norm;

    J(1,1) = dN1_dd1(d(1),d(2));
    J(1,2) = dN1_dd2(d(1),d(2));
    J(2,1) = dN2_dd1(d(1),d(2));
    J(2,2) = dN2_dd2(d(1),d(2));

    while true
        Deltad = -J\R;

        G0 = Deltad' * R;

        d = d + Deltad;
        iter_step = iter_step + 1;

        if iter_step == 15
            break;
        end

        R(1) = R1(d(1),d(2),n);
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));

        G = Deltad' * R;

        s = 1;
        % Line Search
        % find a s such that abs(G) <= stol * abs(G0)
        if abs(G) > stol * abs(G0)
            linmax = 10;
            smax = 16;
            sb = 0;
            sa = s;
            Gb = G0;
            Ga = G;

            % temporary variables
            d_temp = d;
            R_temp = zeros(2,1);

            % find bracket on zero
            % if Ga * Gb > 0, [sb,sa] no zero point
            % then enlarge the bracket to [sa, 2sa]
            % compute the corresponding R_temp and Gb, Ga 
            % in new bracket [sa,2sa]
            while (Ga * Gb > 0 && sa < smax)
                sb = sa;
                sa = 2*sa;
                Gb = Ga;
                d_temp = d_temp + sa * Deltad;
                R_temp(1) = R1(d_temp(1),d_temp(2),n);
                R_temp(2) = R2(d_temp(1),d_temp(2),n);
                Ga = Deltad' * R_temp;
            end

            step = sa; 
            G = Ga;
            lin_step = 0;
            % now we have already found the bracket [sb,sa]
            % which contains the zero point. Ga * Gb < 0

            % Illinois algorithm to find zero
            % while still abs(G) > stol * abs(G0), criteria is not satisfied
            while (lin_step <= linmax && Ga * Gb < 0 && ...
                    abs(G) > stol * abs(G0) ) %|| abs(sb-sa) > stol * 0.5 * (sb+sa)))
                step = sa - Ga * (sa - sb) / (Ga - Gb);   % linear interpolation
                d_temp = d_temp + step * Deltad;
                R_temp(1) = R1(d_temp(1),d_temp(2),n);
                R_temp(2) = R2(d_temp(1),d_temp(2),n);
                G = Deltad' * R_temp;

                if G * Ga > 0
                    Gb = 0.5 * Gb;
                else
                    sb = sa;
                    Gb = Ga;
                end
                sa = step;
                Ga = G;
                lin_step = lin_step + 1;
            end
            s = step;
            % update the new updates
            d = d - Deltad;
            d = d + s * Deltad;
    
            R(1) = R1(d(1),d(2),n);
            R(2) = R2(d(1),d(2),n);
            R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
        end

        if R_norm <= tol * R_norm0
            break;
        end
    end  % end of iteration
    iter_num(:,n+1) = iter_step;
    d_sol(:,n+1) = d;
end % end of load step


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

    R(1) = R1(d(1),d(2),n);
    R(2) = R2(d(1),d(2),n);
    R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
    R_norm0 = R_norm;

    J(1,1) = dN1_dd1(d(1),d(2));
    J(1,2) = dN1_dd2(d(1),d(2));
    J(2,1) = dN2_dd1(d(1),d(2));
    J(2,2) = dN2_dd2(d(1),d(2));

    while true
        Deltad = -J\R;

        G0 = Deltad' * R;

        d = d + Deltad;
        iter_step = iter_step + 1;

        if iter_step == 15
            break;
        end

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
    end  % end of iteration
    iter_num(:,n+1) = iter_step;
    d_sol(:,n+1) = d;
end % end of load step


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



