% clear the memory and the screen
clear; clc;
% =========================================================================
% This is a numerical case for snap-through buckling of a hinged
% right-angle frame using Euler-Bernoulli beam element
%         |
%         V
%   ----------------------------o
%   |                          / \
%   |                         o---o
%   |
%   |
%   |
%   |
%   |
%   |
%   o
%  / \
% o---o
% Consider axial deformation, 
% element nodal displacement column: d^e = [u1 w1 theta1 u2 w2 theta2]^T
% element nodal force column: f^e = [fx1 fz1 m1 fx2 fz2 m2]^T
% element stiffness matrix: k^e : 6 x 6


% =========================================================================
% Problem definition
% exact solution
% exact   = @(x) (x - 1).^5;
% exact_x = @(x) 5 * (x - 1).^4;
% exact_xx = @(x) 20 * (x - 1).^3;
% exact_xxx = @(x) 60 * (x - 1).^2;
% exact_xxxx = @(x) 120 * (x - 1);

% model data
L = 120;
E = 7.2e6;
I = 3.0;
A = 6;
nu = 0.3;

f = @(x) 0;
F_applied = -40000
;
M = 0;
Q = 0;
g1 = 0;
g2 = 0;
g3 = 0;
g4 = 0;
% =========================================================================

% parameters of the FEM
n_sd = 2;              % space dimension
n_el = 10;             % number of elements
n_en = 2;              % number of element nodes
n_ed = 3;              % number of element degrees of freedom (per node)
deg  = n_en - 1;       % polynomial degree
n_np = n_el * deg + 1; % number of points
n_eq = n_np * n_ed - 4;       % number of equations
n_ee = n_ed * n_en;             % number of equations of an element

% quadrature rule
n_int = 2;              % number of quadrature points

% =========================================================================
% Generate the mesh
% generate the coordinates of the nodal points
x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);

hh     = L / (n_el / 2);
x_coor_ver = 0;
x_coor_hor = 0: hh/deg : L;
y_coor_ver = 0: hh/deg : L;
y_coor_hor = L;

mid_pt = floor( n_np / 2 ) + 1;    % 刚性转角点

for n = 1:n_np
    if n <= mid_pt
        x_coor(n) = x_coor_ver;
        y_coor(n) = y_coor_ver(n);
    else
        x_coor(n) = x_coor_hor(n - mid_pt + 1);
        y_coor(n) = y_coor_hor;
    end
end


% ID and LM arrays are generated based on the BC info
ID = zeros(n_ed, n_np);
index = 0;
for AA = 1 : n_np
    for ii = 1 : n_ed      
        if (AA == 1 && ii == 1) || (AA == 1 && ii == 2) ... 
        || (AA == n_np && ii == 1) || (AA == n_np && ii == 2)  % || (AA == mid_pt && ii == 3) 
            ID(ii, AA) = 0;
        else
            index = index + 1;
            ID(ii, AA) = index;
        end
    end
end

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

% Arc length
point_num = 40;     % point numbers
delta_a = 0.3;      % arc-length increment

% arc-length function
% deltad: the total change over the step to current iterate

arclength_fun = @(deltad, deltalambda) deltad' *  deltad + deltalambda * deltalambda - delta_a * delta_a;

% =========================================================================
% generate the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);
% =========================================================================

K = spalloc(n_eq, n_eq, 7*n_eq); % allocate the global stiffness matrix
F = zeros(n_eq, 1);              % allocate the global load vector

x_ele = [0, hh; hh, 2*hh; 2*hh, 3*hh; 3*hh, 4*hh; 4*hh, 5*hh;
         0, hh; hh, 2*hh; 2*hh, 3*hh; 3*hh, 4*hh; 4*hh, 5*hh];

% Assembly of K and F
for ee = 1 : n_el
    k_e = zeros(n_ee, n_ee);         % dimension of element stiffness is n_ee x n_ee
    f_e = zeros(n_ee, 1);            % element force nodes

    % there is only one variable field in each element
    % x_ele = zeros(n_en, 1);
    % 
    % for aa = 1:n_en
    %     AA = IEN(aa, ee);
    %     if AA <= mid_pt
    %         x_ele(aa) = y_coor(AA);
    %     else
    %         x_ele(aa) = x_coor(AA);
    %     end
    % end

    % loop over quadrature points
    for ll = 1 : n_int

        dx_dxi = hh / 2;
        x_l = (hh * xi(ll) + x_ele(ee, 1) + x_ele(ee, 2)) / 2;

        dxi_dx = 1.0 / dx_dxi;
    
        for aa = 1 : n_en
            for ii = 1: n_ed
                pp = n_ed * (aa - 1) + ii;
                if ii == 1
                    f_e(pp) = f_e(pp) + weight(ll) * PolyShape(deg, aa, xi(ll), 0) * f(x_l) * dx_dxi;
                else
                    pp_ind = (n_ed - 1) * (aa - 1) + ii - 1;
                    f_e(pp) = f_e(pp) + weight(ll) * HermiteShape(pp_ind, xi(ll), 0, hh) * f(x_l) * dx_dxi;
                end
                
                for bb = 1 : n_en
                    for jj = 1 : n_ed
                        qq = n_ed * (bb - 1) + jj;
                        if ii == 1 && jj == 1
                            k_e(pp,qq) = k_e(pp,qq) + E * A * weight(ll) * PolyShape(deg, aa, xi(ll), 1) * PolyShape(deg, bb, xi(ll), 1) * dxi_dx;
                        end
                        if ii ~= 1 && jj ~= 1
                            pp_ind = (n_ed - 1) * (aa - 1) + ii - 1;
                            qq_ind = (n_ed - 1) * (bb - 1) + jj - 1;
                            k_e(pp,qq) = k_e(pp,qq) + E * I * weight(ll) * HermiteShape(pp_ind, xi(ll), 2, hh)...
                            * HermiteShape(qq_ind, xi(ll), 2, hh) * dxi_dx * dxi_dx * dxi_dx;
                        end
                    end
                end
            end
        end
    end

    % Compute the rotation matrix
    node1 = [ x_coor(IEN(1, ee)) ,   y_coor(IEN(1, ee))    ];
    node2 = [ x_coor(IEN(n_en, ee)), y_coor(IEN(n_en, ee)) ];
    
    % 计算单元方向向量
    dx = node2(1) - node1(1);
    dy = node2(2) - node1(2);
    L = sqrt(dx^2 + dy^2);
    
    % 计算方向余弦
    cos_theta = dx / L;
    sin_theta = dy / L;
    
    % 构造旋转矩阵
    T = [cos_theta,  sin_theta, 0, 0, 0, 0;
         -sin_theta, cos_theta, 0, 0, 0, 0;
         0,          0,         1, 0, 0, 0;
         0,          0,         0, cos_theta, sin_theta, 0;
         0,          0,         0, -sin_theta, cos_theta, 0;
         0,          0,         0, 0, 0, 1];

    % Trandform the stiffness matrix in terms of the global coordinates
    k_e = T' * k_e * T;


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

    for aa = 1:n_en
        AA = IEN(aa, ee);
        if AA == floor(n_np / 2) + 2
            F(LM(2, AA)) = F_applied;
        end
    end
end

% =========================================================================
% Now we have K and F assembled and we solve the linear system Kd = F
d_temp = K \ F;
detK = det(K);
% Solve linear system using gmres x = gmres(A,b,restart,tol)
% restart = 10000;
% tol = 10e-6;
% d_temp = gmres(K,F,restart,tol);
% =========================================================================


% Generate the full solution vector by inserting back the Dirichlet value
disp = [g1; g2; d_temp(1:end-1); g3; g4; d_temp(end)];
% Extract the displacement component
u_x = zeros(n_np, 1);
u_y = zeros(n_np, 1);
% Extract the angle component
du_dx = zeros(n_np, 1);
ind1 = 0;
ind2 = 0;
ind3 = 0;
for ii = 1 : n_np * n_ed
    if mod(ii, n_ed) == 0
        ind1 = ind1 + 1;
        du_dx(ind1) = disp(ii);
    elseif mod(ii, n_ed) == 1
        ind2 = ind2 + 1;
        u_x(ind2) = disp(ii);
    elseif mod(ii, n_ed) == 2
        ind3 = ind3 + 1;
        u_y(ind3) = disp(ii);
    end
end



% plot the displacement
figure;
plot(x_coor, y_coor,'b-', 'LineWidth', 2);
% xlim([-10,130]);
% ylim([0,130]);
hold on;
% plot displacement
plot(x_coor + u_x, y_coor + u_y,'r-', 'LineWidth', 2);
xlabel("u_x");
ylabel("u_y");
legend('undeformed', 'deformed','Location', 'Best', 'FontSize', 14, 'Box', 'on');
