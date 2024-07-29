% clear the memory and the screen
clear; clc;

% =========================================================================
% Problem definition
% exact solution
exact   = @(x) x.^5;
exact_x = @(x) 5 * x.^4;

f = @(x) -20.0 * x.^3;
g = 1;
h = 0;
% =========================================================================

% parameters of the FEM
n_el  = 5;             % number of elements
n_en  = 2;              % number of element nodes
deg   = n_en - 1;       % polynomial degree
n_np  = n_el * deg + 1; % number of points
n_eq  = n_np - 1;       % number of equations
n_int = 3;              % number of quadrature points

% =========================================================================
% Generate the mesh
% nodal coordinates
hh     = 1 / n_el;
x_coor = 0 : hh/deg : 1;

% IEN
IEN = zeros(n_en, n_el);

for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(aa,ee) = (ee-1) * deg + aa;
  end
end
% =========================================================================

% ID and LM arrays are generated based on the BC info
ID = 1 : n_np;
ID(end) = 0; % Modify ID according to the Dirichlet BC info

LM = ID(IEN);

% generate the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

K = spalloc(n_eq, n_eq, 3*n_eq); % allocate the global stiffness matrix
F = zeros(n_eq, 1);    % allocate the global load vector

% Assembly of K and F
for ee = 1 : n_el
    k_e = zeros(n_en, n_en);
    f_e = zeros(n_en, 1);
    
    x_ele = x_coor(IEN(1:n_en,ee)); % A = IEN(a,e) and x_ele(a) = x_coor(A)

    for ll = 1 : n_int
        dx_dxi = 0.0;
        x_l = 0.0;
        % for aa = 1 : n_en
        %     dx_dxi = dx_dxi + x_ele(aa) * PolyShape(deg, aa, xi(ll), 1);
        %     x_l = x_l + x_ele(aa) * PolyShape(deg, aa, xi(ll), 0);
        % end
        dx_dxi = hh / 2;
        x_l = (hh * xi(ll) + x_ele(1) + x_ele(2)) / 2;
        
        dxi_dx = 1.0 / dx_dxi;
    
        for aa = 1 : n_en
            f_e(aa) = f_e(aa) + weight(ll) * PolyShape(deg, aa, xi(ll), 0) * f(x_l) * dx_dxi;
            for bb = 1 : n_en
                k_e(aa,bb) = k_e(aa,bb) + weight(ll) * PolyShape(deg, aa, xi(ll), 1) * PolyShape(deg, bb, xi(ll), 1) * dxi_dx;
            end
        end
    end

    % Now we need to put element k and f into global K and F
    for aa = 1 : n_en
        PP = LM(aa,ee);
        if PP > 0
            F(PP) = F(PP) + f_e(aa);
            for bb = 1 : n_en
                QQ = LM(bb,ee);
                if QQ > 0
                    K(PP,QQ) = K(PP,QQ) + k_e(aa,bb);
                else
                    F(PP) = F(PP) - k_e(aa,bb) * g;
                end
            end
        end
    end

    if ee == 1
        F(ID(IEN(1,ee))) = F(ID(IEN(1,ee))) + h;
    end
end

% =========================================================================
% Now we have K and F assembled and we solve the linear system Kd = F
d_temp = K \ F;

% Solve linear system using gmres x = gmres(A,b,restart,tol)
% restart = 10000;
% tol = 10e-6;
% d_temp = gmres(K,F,restart,tol);
% =========================================================================


% Generate the full solution vector by inserting back the Dirichlet value
disp = [d_temp; g];

% plot the solution
figure;
X_h = 0: hh/deg :1;
Y_h = disp;
% subplot(2,1,1)
plot(X_h, Y_h,'b-', 'LineWidth', 2);
hold on;
% subplot(2,1,2)
X = 0:0.01:1;
Y = exact(X);
plot(X, Y,'r-', 'LineWidth', 2);
xlabel("X");
ylabel("Displacement");
legend('有限元解', '精确解','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% calculate the error
nqp = 10; % we need more points 
[xi, weight] = Gauss(nqp, -1, 1);

L2_top = 0.0; L2_bot = 0.0; H1_top = 0.0; H1_bot = 0.0;

for ee = 1 : n_el
    for ll = 1 : nqp
        x_ele = x_coor( IEN(1:n_en, ee) );
        u_ele = disp(   IEN(1:n_en, ee) );

        x_l = 0.0; uh = 0.0; dx_dxi = 0.0; du_dxi = 0.0;
        for aa = 1 : n_en
            x_l    = x_l    + x_ele(aa) * PolyShape(deg, aa, xi(ll), 0);
            uh     = uh     + u_ele(aa) * PolyShape(deg, aa, xi(ll), 0);
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(deg, aa, xi(ll), 1);
            du_dxi = du_dxi + u_ele(aa) * PolyShape(deg, aa, xi(ll), 1);
        end

        dxi_dx = 1.0 / dx_dxi;

        L2_top = L2_top + weight(ll) * ( uh - exact(x_l) )^2 * dx_dxi;
        L2_bot = L2_bot + weight(ll) * exact(x_l)^2 * dx_dxi;

        H1_top = H1_top + weight(ll) * (du_dxi*dxi_dx - exact_x(x_l))^2 * dx_dxi;
        H1_bot = H1_bot + weight(ll) * exact_x(x_l)^2 * dx_dxi;
    end
end

L2_top = sqrt(L2_top); L2_bot = sqrt(L2_bot);

H1_top = sqrt(H1_top); H1_bot = sqrt(H1_bot);

L2_error = L2_top / L2_bot;
H1_error = H1_top / H1_bot;




% eof