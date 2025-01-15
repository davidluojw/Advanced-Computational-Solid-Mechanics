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
d = zeros(n_eq,1);      % d(i)
uh = [ d; g(omega_r) ];

% line search
stol = 0.5;

counter = 0;
nmax    = 100;
error   = 1.0;

uh_set = zeros(n_np,nmax + 1);

K = AssemblyK(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa);
% F(i)
F = AssemblyF(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);

while counter < nmax && error > 1.0e-8

    % Solve the stiffness matrix problem
    Deltad = K \ F;      % Deltad(i)
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
    
    % F(i+1)
    F = AssemblyF(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);

    error = norm(F);
    counter = counter + 1;
    uh_set(:, counter) = uh;
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