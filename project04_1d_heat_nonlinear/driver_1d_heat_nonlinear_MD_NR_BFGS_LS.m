% This is the nonlinear solver for heat equation with material nonlinearity
% clean the memory and the screen
clear all; clc;

% domain
omega_l = 0.0;
omega_r = 10.0;

% -------------------------------------------------------------------------
% material properties and input data
%kappa = 1.0;

% exact solution
exact = @(x) sin(x);
exact_x = @(x) cos(x);

% nonlinear kappa
fun_kappa  = @(u) 1 + u*u;
fun_dkappa = @(u) 2*u;

f = @(x) -2.0*cos(x)*cos(x)*sin(x) + sin(x)*(sin(x)*sin(x) + 1.0);
h = @(x) -1.0 * fun_kappa( exact(omega_l) )* exact_x(omega_l);
g = @(x) exact(omega_r);
% -------------------------------------------------------------------------

% interpolation degree
pp = 2;

% number of elements
nElem = 10;

% quadrature rule
nqp = 2; %pp + 1;
[qp, wq] = Gauss(nqp, -1, 1);

n_np = nElem * pp + 1; % number of nodal points
n_en = pp + 1;         % number of element nodes

IEN = zeros(n_en, nElem);

for ee = 1 : nElem
    for aa = 1 : n_en
        IEN(aa, ee) = (ee - 1) * pp + aa;
    end
end

% mesh is assumbed to have uniform size hh
hh = (omega_r - omega_l) / nElem;

x_coor = omega_l : (hh/pp) : omega_r;

% setup ID array based on the boundary condition
ID = 1 : n_np;
ID(end) = 0;

LM = ID (IEN);

% Setup the stiffness matrix and load vector
% number of equations equals the number of nodes minus the number of
% Dirichlet nodes
n_eq = n_np - 1;

% initial guess
d = zeros(n_eq,1);
uh = [ d; g(omega_r) ];

% line search
stol = 0.5;

counter = 0;
nmax    = 100;
error   = 1.0;

% update vectors: V and W
V = zeros(n_eq, nmax + 1);
W = zeros(n_eq, nmax + 1);

uh_set = zeros(n_np,nmax + 1);
error_set = zeros(nmax+1, 1);
s_set = zeros(nmax+1,1);

counter = counter + 1;

% get the difference between first two residual
% Delta F = F(i+1) - F(i)

F = AssemblyF(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);
error = norm(F);

F_old = F;   % F(i)
d_old = d;   % d(i)

K = AssemblyK(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa);
uh_set(:, counter) = uh;
error_set(counter,:) = error;

Deltad = K \ F;      % Deltad(i)

%Line Search ---------------------------------
d_temp = d;          % d(i)
G0 = Deltad' * F;    % G0 = Deltad(i) * F(i)

d = d + Deltad;  % d(i+1) = d(i) + Deltad(i)
uh = [ d ; g(omega_r) ];

% F(i+1)
F = AssemblyF(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);

G = Deltad' * F;   % G = Deltad(i) * F(i+1)

% d(i), DEltad(i)
s = LineSearch_heat(G0,G,Deltad,d_temp,stol,omega_r,g,pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);

% update the new updates
d = d_temp + s * Deltad;   % d(i+1) = d(i) + s*Deltad(i)
uh = [d; g(omega_r) ];

Deltad = d - d_temp;       % Delta(i) = d(i+1) - d(i)  

%Line Search ---------------------------------

% d = d + Deltad;      %d(i+1) = d(i) + Deltad(i)
% uh = [ d; g(omega_r) ];

counter = counter + 1;

% This is the new residual F(i+1)
F = AssemblyF(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);
error = norm(F);

uh_set(:, counter) = uh;
error_set(counter,:) = error;
s_set(counter,:) = s;

% BFGS start, solve K_bar * Delta d(i+1) = F(i+1)
% where Delta d(i+1) = d(i+2) - d(i+1)
% actually, BFGS is to find d(i+2)
while counter < nmax
    F_tilde = F;    % F_tilde(0) = F(i+1)

    % compute update vectors V and W
    % V(i-k+1) = Deltad(i-k+1) / (Deltad(i-k+1) * DeltaF(i-k+1))
    % DeltaF(i-k+1) = F(i-k+2) - F(i-k+1)
    % Deltad(i-k+1) = d(i-k+2) - d(i-k+1)

    % k = 1, i-k+1 = i, i-k+2 = i+1
    DeltaF_k = F_tilde - F_old;     % Delta F(i-k+1) = Delta F(i) = F(i+1) - F(i)
    Deltad_k = Deltad;              % Delta d(i-k+1) = Delta d(i) 

    % V(i-k+1) = Deltad(i-k+1) / (Deltad(i-k+1) * DeltaF(i-k+1))
    V(:,counter) = Deltad_k / (Deltad_k' * DeltaF_k);

    % alpha = (-Delta F(i-k+1) * Delta d(i-k+1) / (F(i-k+1) * Delta d(i-k+1)) )^(1/2)
    %       = (-Delta F(i) * Delta d(i) / (F(i) * Delta d(i) ) )^(1/2)
    alpha = sqrt(  abs( - s * DeltaF_k' * Deltad_k / (F_old' * Deltad_k)  ) ) ;

    % W(i-k+1) = alpha * F(i-k+1) - Delta F(i-k+1)
    %          = alpha * F(i) - Delta F(i)
    W(:,counter) = -alpha * F_old - DeltaF_k;

    % Right-side updates
    for k = 1:counter
        F_tilde = F_tilde + (V(:, counter - k + 1)'* F_tilde) * W(:,counter - k + 1);
    end

    % Solve intermediate system
    % Deltad_tilde(0) = K^(-1) * F_tilde(i)
    Deltad_tilde = K \ F_tilde;

    % Left-side update
    for k = 1:counter
        Deltad_tilde = Deltad_tilde + (W(:, k)' * Deltad_tilde) * V(:, k);
    end
    d_old = d;             % d_old = d(i+1)
    F_old = F;             % F_old = F(i+1)

    %Line Search ---------------------------------
    d_temp = d;          % d(i+1)
    G0 = Deltad_tilde' * F;    % G0 = Deltad_tilde(i) * F(i+1) = Deltad(i+1)*F(i+1)

    d = d + Deltad_tilde;  % d(i+2) = d(i+1) + Deltad(i+1)
    uh = [ d ; g(omega_r) ];

    % F(i+2)
    F = AssemblyF(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);

    G = Deltad_tilde' * F;   % G = Deltad(i+1) * F(i+1)

    % d(i+1), Deltad(i+1)
    s = LineSearch_heat(G0,G,Deltad_tilde,d_temp,stol,omega_r,g,pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);

    % update the new updates
    d = d_temp + s * Deltad_tilde;   % d(i+2) = d(i+1) + s*Deltad(i+1)
    uh = [d; g(omega_r) ];

    Deltad = d - d_temp;       % Delta(i+1) = d(i+2) - d(i+1)  
    %Line Search ---------------------------------

    % d = d + Deltad_tilde;  % d(i+2) = d(i+1) + Deltad_tilde
    % uh = [ d ; g(omega_r) ];
    % Deltad = d - d_old;    % Deltad(i+1) = d(i+2) - d(i+1)

    counter = counter + 1;

    % F(i+2)
    F = AssemblyF(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);
    error = norm(F);

    if error <= 1.0e-8
        uh_set(:, counter) = uh;
        error_set(counter, :) = error;
        s_set(counter,:) = s;
        break;
    end

    uh_set(:, counter) = uh;
    error_set(counter,:) = error;
    s_set(counter,:) = s;
end

% plot the solution
figure;
X_h = omega_l: hh/pp :omega_r;
for i=1:counter
    Y_h = uh_set(:,i);
    h_fem = plot(X_h, Y_h,'b-', 'LineWidth', 2);
    hold on;
end
X = omega_l:0.01:omega_r;
Y = exact(X);
h_exact = plot(X, Y,'r--', 'LineWidth', 2);

legend([h_fem, h_exact], {'FEM', 'EXACT'}, 'Location', 'Best', 'FontSize', 14, 'Box', 'on');
xlabel("X");
ylabel("Temperature");

figure;
X_h = omega_l: hh/pp :omega_r;
Y_h = uh;
plot(X_h, Y_h,'b-', 'LineWidth', 2);
% subplot(2,1,1)
xlabel("X");
ylabel("Temperature");
hold on;
% subplot(2,1,2)
X = omega_l:0.01:omega_r;
Y = exact(X);
plot(X, Y,'r--', 'LineWidth', 2);
xlabel("X");
ylabel("Temperature");
legend('FEM','EXACT','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% Now we do the postprocessing
nqp = 6;
[qp, wq] = Gauss(nqp, -1, 1);

top = 0.0; bot = 0.0;
top_u = 0.0; bot_u = 0.0;
for ee = 1 : nElem
    for qua = 1 : nqp
        x_ele = zeros(n_en, 1);
        u_ele = zeros(n_en, 1);
        for aa = 1 : n_en
            x_ele(aa) = x_coor(IEN(aa, ee));
            u_ele(aa) = uh(IEN(aa, ee));
        end
        
        x = 0.0; dx_dxi = 0.0; u = 0.0; duh_dxi = 0.0;
        for aa = 1 : n_en
            x = x + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            u = u + u_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            dx_dxi = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
            duh_dxi = duh_dxi + u_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
        end

        top_u = top_u + wq(qua) * ( u - exact(x) )^2 * dx_dxi;
        bot_u = bot_u + wq(qua) * exact(x)^2 * dx_dxi;
        
        dxi_dx = 1.0 / dx_dxi;
        
        top = top + wq(qua) * (duh_dxi * dxi_dx - exact_x(x))^2 * dx_dxi;
        bot = bot + wq(qua) * exact_x(x)^2 * dx_dxi;
    end
end

top = sqrt(top);
bot = sqrt(bot);

% error = top / bot;

top_u = sqrt(top_u);
bot_u = sqrt(bot_u);
error_u = top_u / bot_u;

% EOF