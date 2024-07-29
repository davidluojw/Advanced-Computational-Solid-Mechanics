% clear the memory and the screen
clear; clc;

% =========================================================================
% Problem definition
% exact solution
exact   = @(x) (x - 1).^5;
exact_x = @(x) 5 * (x - 1).^4;
exact_xx = @(x) 20 * (x - 1).^3;
exact_xxx = @(x) 60 * (x - 1).^2;
exact_xxxx = @(x) 120 * (x - 1);

E = 1.0;
I = 1.0;

f = @(x) 120.0 * E * I * (x - 1);
M = -20 * E * I;
Q = 60 * E * I;
g1 = 0;
g2 = 0;
% =========================================================================

% parameters of the FEM
n_el  = 4;             % number of elements
n_en  = 2;              % number of element nodes
n_ed = 2;               % number of element degrees of freedom (per node)
n_np  = n_el * (n_en - 1) + 1; % number of points
n_eq  = n_np * n_ed - 2;     % number of equations
n_int = 3;              % number of quadrature points
n_ee = n_ed * n_en;     % number of element equations

% =========================================================================
% Generate the mesh
% nodal coordinates
hh     = 1 / n_el;
x_coor = 0 : hh/(n_en - 1) : 1;

% ID and LM arrays are generated based on the BC info
ID = zeros(n_ed, n_np);
for AA = 1 : n_np
    for ii = 1 : n_ed
        ID(ii, AA) = (AA - 1) * n_ed + ii;
    end
end
% Modify ID according to the Dirichlet BC info
ID(n_ed - 1, n_np) = 0; 
ID(n_ed, n_np) = 0; 


% IEN
IEN = zeros(n_en, n_el);
for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(aa,ee) = (ee-1) * (n_en - 1) + aa;
  end
end

%LM
LM = zeros(n_ee, n_el);
for ee = 1 : n_el
    for aa = 1 : n_en
        for ii = 1 : n_ed
            pp = n_ed * (aa - 1) + ii;
            LM(pp, ee) = ID(ii, IEN(aa, ee));
        end
    end
end


% =========================================================================
% generate the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);
% =========================================================================

K = spalloc(n_eq, n_eq, 7*n_eq); % allocate the global stiffness matrix
F = zeros(n_eq, 1);    % allocate the global load vector

% Assembly of K and F
for ee = 1 : n_el
    k_e = zeros(n_ee, n_ee);     % dimension of element stiffness is n_ee x n_ee
    f_e = zeros(n_ee, 1);
    
    x_ele = x_coor(IEN(1:n_en,ee)); % A = IEN(a,e) and x_ele(a) = x_coor(A)

    for ll = 1 : n_int

        % for aa = 1 : n_en
        %     dx_dxi = dx_dxi + x_ele(aa) * PolyShape(deg, aa, xi(ll), 1);
        %     x_l = x_l + x_ele(aa) * PolyShape(deg, aa, xi(ll), 0);
        % end
        dx_dxi = hh / 2;
        x_l = (hh * xi(ll) + x_ele(1) + x_ele(2)) / 2;

        dxi_dx = 1.0 / dx_dxi;
    
        for aa = 1 : n_en
            for ii = 1: n_ed
                pp = n_ed * (aa - 1) + ii;
                f_e(pp) = f_e(pp) + weight(ll) * HermiteShape(pp, xi(ll), 0, hh) * f(x_l) * dx_dxi;
                for bb = 1 : n_en
                    for jj = 1 : n_ed
                        qq = n_ed * (bb - 1) + jj;
                        k_e(pp,qq) = k_e(pp,qq) + weight(ll) * HermiteShape(pp, xi(ll), 2, hh)...
                            * HermiteShape(qq, xi(ll), 2, hh) * dxi_dx * dxi_dx * dxi_dx;
                    end
                end
            end
        end
    end

    % Now we need to put element k and f into global K and F
    for aa = 1 : n_en
        for ii = 1 : n_ed
            pp = n_ed * (aa - 1) + ii;
            PP = LM(pp,ee);
            if PP > 0
                F(PP) = F(PP) + f_e(pp);
                for bb = 1 : n_en
                    for jj = 1 : n_ed
                        qq = n_ed * (bb - 1) + jj;
                        QQ = LM(qq,ee);
                        if QQ > 0
                            K(PP,QQ) = K(PP,QQ) + k_e(pp,qq);
                        else
                            F(PP) = F(PP) - k_e(pp,qq) * g1 - k_e(pp,qq) * g2;
                        end
                    end
                end
            end
        end
    end

    if ee == 1
        F(ID(1,IEN(1,ee))) = F(ID(1, IEN(1,ee))) + Q;
        F(ID(2,IEN(1,ee))) = F(ID(2, IEN(1,ee))) - M;
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
disp = [d_temp; g1; g2];
% Extract the displacement component
u_h = zeros(n_np * n_ed / 2, 1);
% Extract the angle component
du_dx = zeros(n_np * n_ed / 2, 1);
for ii = 1 : n_np * n_ed
    if mod(ii, 2) == 1
        u_h(floor(ii / 2) + 1) = disp(ii);
    else
        du_dx(floor(ii / 2)) = disp(ii);
    end
end



% plot the solution
figure;
X_h = 0: hh / (n_en - 1) :1;
Y_h = u_h;
plot(X_h, Y_h,'b-', 'LineWidth', 2);
hold on;
X = 0:0.01:1;
Y = exact(X);
plot(X, Y,'r-', 'LineWidth', 2);
xlabel("X");
ylabel("Displacement");
legend('有限元解-位移', '精确解-位移','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% plot the angle
figure;
X_h = 0: hh / (n_en - 1) :1;
Y_h = du_dx;
plot(X_h, Y_h,'b-', 'LineWidth', 2);
hold on;
X = 0:0.01:1;
Y = exact_x(X);
plot(X, Y,'r-', 'LineWidth', 2);
xlabel("X");
ylabel("Angle");
legend('有限元解-转角', '精确解-转角','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% calculate the error
nqp = 10; % we need more points 
[xi, weight] = Gauss(nqp, -1, 1);

L2_top = 0.0; L2_bot = 0.0; H1_top = 0.0; H1_bot = 0.0;

for ee = 1 : n_el
    for ll = 1 : nqp
        x_ele = x_coor( IEN(1:n_en, ee) );

        uh = 0.0;            % displacement at xi(ll)
        duh_dxi = 0.0;        % derivative of displacement at xi(ll)

        dx_dxi = hh / 2;
        x_l = (hh * xi(ll) + x_ele(1) + x_ele(2)) / 2;

        if ee == n_el                         % the last element, which do not invole the g-nodes
            d_ele = disp( LM(1:n_ed, ee)  );
            d_x_ele = disp( LM(1:n_ed, ee) );
            for ii = 1 : n_ed
                uh = uh + d_ele(ii) * HermiteShape(ii, xi(ll), 0, hh);
                duh_dxi = duh_dxi + d_ele(ii) * HermiteShape(ii, xi(ll), 1, hh);
            end

        else
            d_ele = disp( LM(1:n_en * n_ed, ee)  );
            d_x_ele = disp( LM(1:n_en * n_ed, ee) );
            for aa = 1 : n_en
                for ii = 1 : n_ed
                    pp = n_ed * (aa - 1) + ii;
                    uh = uh + d_ele(pp) * HermiteShape(pp, xi(ll), 0, hh);
                    duh_dxi = duh_dxi + d_ele(pp) * HermiteShape(pp, xi(ll), 1, hh);
                end
            end
        end
        

        

        dxi_dx = 1.0 / dx_dxi;

        L2_top = L2_top + weight(ll) * ( uh - exact(x_l) )^2 * dx_dxi;
        L2_bot = L2_bot + weight(ll) * exact(x_l)^2 * dx_dxi;

        H1_top = H1_top + weight(ll) * (duh_dxi*dxi_dx - exact_x(x_l))^2 * dx_dxi;
        H1_bot = H1_bot + weight(ll) * exact_x(x_l)^2 * dx_dxi;
    end
end

L2_top = sqrt(L2_top); L2_bot = sqrt(L2_bot);

H1_top = sqrt(H1_top); H1_bot = sqrt(H1_bot);

L2_error = L2_top / L2_bot;
H1_error = H1_top / H1_bot;




% eof