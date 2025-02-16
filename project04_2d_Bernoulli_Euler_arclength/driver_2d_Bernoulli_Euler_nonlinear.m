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

% model data
L = 120;
E = 7.2e6;
I = 2;
A = 6;
nu = 0.3;

f = @(x) 0;
F_applied = -27000;
P_bar = 0;
M_bar = 0;
Q_bar = 0;
n_gBC = 4;
g_BC = zeros(n_gBC, 1);
% =========================================================================

% parameters of the FEM
n_sd = 2;              % space dimension
n_el = 20;             % number of elements
n_en = 2;              % number of element nodes
n_ed = 3;              % number of element degrees of freedom (per node)
deg  = n_en - 1;       % polynomial degree
n_np = n_el * deg + 1; % number of points
n_eq = n_np * n_ed - n_gBC;       % number of equations
n_ee = n_ed * n_en;           % number of equations of an element

% quadrature rule
n_int = 2;              % number of quadrature points

% =========================================================================

% Model
Model = struct( 'n_sd', n_sd, 'n_el', n_el, 'n_en', n_en, 'n_ed', n_ed, 'deg', deg, 'n_np', n_np, ...
                'n_eq', n_eq, 'n_ee', n_ee, 'n_int', n_int, 'E', E, 'I', I, 'A', A);

% =========================================================================
% Generate the mesh
% generate the coordinates of the nodal points
x_coor     = zeros(n_np, 1);   % Absolute x coordinate in global
y_coor     = zeros(n_np, 1);   % Absolute y coordinate in global
theta_coor = zeros(n_np, 1);   % Absolute theta coordinate in global

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
% =========================================================================

% Geometry
Geometry = struct('x_coor', x_coor, 'y_coor', y_coor, 'theta_coor', theta_coor, ...
                   'mid_pt', mid_pt, 'ID', ID, 'IEN', IEN, 'LM', LM, 'L', L);


% =========================================================================

% Boundary conditions
BC = struct('f', f, 'F_applied', F_applied, 'g_BC', g_BC, 'P_bar', P_bar, 'Q_bar', Q_bar, 'M_bar', M_bar);

% =========================================================================
% Arc length
point_num = 50;     % point numbers
delta_a = 0.03;      % arc-length increment

% =========================================================================
% generate the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);
% =========================================================================


% =========================================================================
% Newton-Raphson Algorithm with using arc-length method
% initial point value
disp = zeros(n_np * n_ed,1);
d    = zeros(n_eq, 1);
% the index for d, which is used in the calculation of the matrix
d_ind = [3:(n_np * n_ed-3), n_np * n_ed];
d     = disp(d_ind);
lambda = 0.0;
% tolerance
tol = 1.0e-8;

% increment of the load step
deltad = zeros(n_eq,1);
deltalambda = 0;

% increment of the iteration step
Deltad = zeros(n_eq,1);
Deltalambda = 0;

% initial increment of lamdba
Deltalambda_init = 0.01;


% =========================================================================
% solver
Solver = struct('point_num', point_num, 'delta_a', delta_a, 'xi', xi, 'weight', weight, ... 
                'Deltalambda_init', Deltalambda_init);

% =========================================================================

Results = struct('d_step', [], 'disp_step',[], 'lambda_step', [], 'd_iter', [], 'disp_iter', [], 'lambda_iter', [], 'd_ind', d_ind);
Results.d_step(:, 1)    = d;
Results.disp_step(:, 1) = disp;
Results.lambda_step(1)  = lambda;
Results.d_iter(:, 1)    = d;
Results.disp_iter(:, 1) = disp;
Results.lambda_iter(1)  = lambda;

% =========================================================================

% F_ext 
F_ext = AssemblyF(Model, Geometry, BC, Solver);
F_ext_norm = norm(F_ext);

Force_ext = [0 ; 0; F_ext(1:end-1); 0; 0; F_ext(end)];

% =========================================================================

Node(n_np, 1) = struct('x', [], 'y', [], 'P', [], 'Q', [], 'M', [],'dof',[]);

for ii = 1: n_np
    Node(ii).x   = x_coor(ii);
    Node(ii).y   = y_coor(ii);
    Node(ii).P   = Force_ext(n_ed * (ii - 1) + 1);
    Node(ii).Q   = Force_ext(n_ed * (ii - 1) + 2);
    Node(ii).M   = Force_ext(n_ed * (ii - 1) + 3);
    Node(ii).dof = [n_ed * (ii - 1) + 1, n_ed * (ii - 1) + 2, n_ed * (ii - 1) + 3 ];
end


% =========================================================================
% Updated Lagrangian 
% Elem: save all the information of the element at the n-th loading step
Elem(n_el,1) = struct('Node1',[], 'Node2', [], 'eleLength',[],'eleAngle',[], ...
                      'eleFint',[], 'eleTangentK',[],'eleElasticK',[]);

% Initialization
for ee = 1:n_el
    % element nodes
    Node1 = Node(IEN(1, ee));
    Node2 = Node(IEN(2, ee));

    % local coordinates, x along axial, y along normal to axial, z normal to the surface
    vl = [Node2.x - Node1.x, Node2.y - Node1.y, 0];
    vx = vl / norm(vl);
    vz = [0, 0, 1];
    vy = cross(vz, vx);

    Elem(ee).Node1              = Node1;
    Elem(ee).Node2              = Node2;
    Elem(ee).eleLength          = norm(vl);

    Elem(ee).eleAngle           = atan2(Node2.y - Node1.y, Node2.x - Node1.x);  % thisi the axial line angle, not the deformation angle rotatio
    
    Elem(ee).eleFint            = zeros(n_ee,1);
    Elem(ee).eleTangentK        = zeros(n_ee, n_ee);
    Elem(ee).eleElasticK        = zeros(n_ee, n_ee);
end

Elem_n = Elem;



% =========================================================================

% arc-length function
% deltad: the total change over the step to current iterate
arclength_fun = @(deltaU, deltaLam) deltaU' *  deltaU + deltaLam * deltaLam * F_ext_norm * F_ext_norm - delta_a * delta_a;

% Arc-length process
for n = 1:point_num+1

    d_n = d;
    lambda_n = lambda;
    disp_n = disp;

    % Stored the last load step configuration
    for ee = 1:n_el
        AA = [IEN(1, ee), IEN(2,ee)]';
        PP = [ n_ed * (ee - 1) + 1, n_ed * (ee - 1) + 2, n_ed * (ee - 1) + 3, ...
               n_ed * ee + 1,       n_ed * ee + 2,       n_ed * ee + 3];
        
        Elem_n(ee).eleLength          = Elem(ee).eleLength;
        Elem_n(ee).eleAngle           = Elem(ee).eleAngle;  
        Elem_n(ee).eleFint            = Elem(ee).eleFint;
        Elem_n(ee).eleTangentK        = Elem(ee).eleTangentK ;
        Elem_n(ee).eleElasticK        = Elem(ee).eleElasticK ;
    end


    % K(d_n+1^(0))
    [K0, Elem] = AssemblyK_SecondOrder(Model, Geometry, Solver, Elem, Elem_n, disp, disp_n);   % updated Elem.eleElasticK and Elem.eleTangentK
    detJ = det(K0);   % det K(d_n+1^(0))

    q_bar = K0 \ F_ext;          % q_bar_n+1^(0) = K(d_n+1^(0))^(-1) * F_ext    

    if n == 1
        % 手动设置
        s = 1;
        % Deltalambda = sqrt( delta_a * delta_a / (F_ext_norm * F_ext_norm  + q_bar' * q_bar));  % Deltalambda^(0)
        Deltalambda = Deltalambda_init;
        n2 = q_bar' * q_bar;                    % qbar_0^(0) * qbar_0^(0): used for GSP calculation
    else
        % Generalized Stiffness Parameter
        GSP = n2 / (q_bar_n' * q_bar);        % GSP = n2 / qbar_n^(0) * qbar_n+1^(0)

        % if GSP is negative, convert the sign of the increment
        if GSP < 0
            s = -s;
        end

        % adopt constant arclength - spherical
        % Deltalambda_n+1^(0) = s * sqrt(delta_a^2 / (Fextnorm^2 + qbar_n+1^(0)^2))
        Deltalambda = s * sqrt( (deltad'*deltad + deltalambda^2 * F_ext_norm * F_ext_norm) / (F_ext_norm * F_ext_norm + q_bar' * q_bar)); 

    end

    Deltad = Deltalambda * q_bar;  % Deltad_n+1^(0) = Deltalambda_n+1^(0)* q_bar_n+1^(0)

    % Update points
    d = d + Deltad;   % d_n+1^(1) = d_n+1^(0) + Deltad_n+1^(0)
    lambda = lambda + Deltalambda;   % lambda_n+1^(1) = lambda_n+1^(0) + Deltalambda^(0)
    disp = [g_BC(1); g_BC(2); d(1:end-1); g_BC(3); g_BC(4); d(end)];

    deltad      = d - d_n;                  % deltad_n+1^(0)
    deltalambda = lambda - lambda_n;    % deltalambda_n+1^(0)
    deltadisp   = disp - disp_n;

    % Stored the current load step initial increment for the next load step
    Deltad_n      = Deltad;
    Deltalambda_n = Deltalambda;
    q_bar_n       = q_bar;    % will be used for the next GSP calculation
    
    % =====================================================================
    % Starting iterative process
    iter_step = 0;   % iteration step
    % N-R iteration with constraint surface
    while true

        % R(d_n+1^(i), lambda_n+1^(i))
        [N, Elem] = AssemblyN_SecondOrder(Model, Geometry, Solver, Elem, Elem_n, disp, disp_n);
        R = lambda * F_ext - N;
        R_norm = norm(R);

        % Store iteration results
        Results.d_iter(:, size(Results.d_iter,2)+1)       = d;
        Results.disp_iter(:, size(Results.disp_iter,2)+1) = disp;
        Results.lambda_iter(size(Results.lambda_iter,2)+1)  = lambda;

        if R_norm <= tol * F_ext_norm
            break;
        end

        iter_step = iter_step + 1;
    
        if iter_step == 100
            break;
        end

        [K0, Elem] = AssemblyK_SecondOrder(Model, Geometry, Solver, Elem, Elem_n, disp, disp_n); 

        q_bar = K0 \ F_ext;     % qbar_n+1^(i) = K(d_n+1^(i))^(-1) * F_ext 

        Deltad_bar = K0 \ R;    % Deltad_bar_n+1^(i) = K(d_n+1^(i))^(-1) * R(d_n+1^(i), lambda_n+1^(i))

        % Constant arc-length spherical
        a = q_bar' * q_bar + F_ext_norm * F_ext_norm;
        b = q_bar' * (Deltad_bar + deltad) + deltalambda * F_ext_norm * F_ext_norm;
        c = Deltad_bar' * (Deltad_bar + 2 * deltad);
        s_iter = sign(deltad' * q_bar);
        Deltalambda = -b/a + s_iter * sqrt((b/a)^2 -  c/ a);  % Deltalambda_n+1^(i)

        %为何会出现复数
        if (~isreal(Deltalambda))
            % conv = -1;
            break;
        end

        % Deltad^(i) = Deltalambda^(i) * q_bar + Deltad_bar
        Deltad = Deltalambda * q_bar + Deltad_bar;

        % Update total values of load ratio and displacements
        d = d + Deltad;                 % d_n+1^(i+1) = d_n+1^(i) + Deltad^(i)
        lambda = lambda + Deltalambda;  % lambda_n+1^(i+1) = lambda_n+1^(i) + Deltalambda^(i)
        disp = [g_BC(1); g_BC(2); d(1:end-1); g_BC(3); g_BC(4); d(end)];

        deltad      = d - d_n;                  % deltad_n+1^(i+1) = d_n+1^(i+1) - d_n
        deltalambda = lambda - lambda_n;        % deltalambda_n+1^(i+1) = lambda_n+1^(i+1) - lambda_n
        deltadisp   = disp - disp_n;           

    end

    iter_set(:,n) = iter_step;

    % Store the solutions: equilibrium configuration
    Results.d_step(:, n+1)    = d;
    Results.disp_step(:, n+1) = disp;
    Results.lambda_step(n+1)  = lambda;

    if (lambda >= 1.0)
        break;
    end
 

end

% point_num = n;


% Draw structure (initial and deformed configurations)
figure;
size_load = size(Results.disp_step,2);
for col=1:size_load
    for ee = 1:n_el
        x1  = Elem(ee).Node1.x + Results.disp_step(Elem(ee).Node1.dof(1),col);
        y1  = Elem(ee).Node1.y + Results.disp_step(Elem(ee).Node1.dof(2),col);
        x2  = Elem(ee).Node2.x + Results.disp_step(Elem(ee).Node2.dof(1),col);
        y2  = Elem(ee).Node2.y + Results.disp_step(Elem(ee).Node2.dof(2),col);
    
        x   = [x1 x2];
        y   = [y1 y2];
    
        scatter(x,y,15,'filled','black');
        hold on
        plot(x,y,'black');
        hold on
    end

    

end
hold off;

figure;
size_iter = size(Results.disp_iter,2);
for col=1:floor(size_iter / 20):size_iter
    for ee = 1:n_el
        x1  = Elem(ee).Node1.x + Results.disp_iter(Elem(ee).Node1.dof(1),col);
        y1  = Elem(ee).Node1.y + Results.disp_iter(Elem(ee).Node1.dof(2),col);
        x2  = Elem(ee).Node2.x + Results.disp_iter(Elem(ee).Node2.dof(1),col);
        y2  = Elem(ee).Node2.y + Results.disp_iter(Elem(ee).Node2.dof(2),col);

        x   = [x1 x2];
        y   = [y1 y2];

        scatter(x,y,15,'filled','black');
        hold on
        plot(x,y,'black');
        hold on
    end
end
hold off;

figure;
for ii = 1:n_np
    if Node(ii).x == 24
        plot(Results.disp_iter(Node(ii).dof(1),:),Results.lambda_iter,'Marker','x','MarkerSize',3);
    end
end

figure;
for ii = 1:n_np
    if Node(ii).x == 24
        plot(Results.disp_iter(Node(ii).dof(2),:),Results.lambda_iter,'Marker','x','MarkerSize',3);
    end
end




% plot3(Result.U_step(Node(nid(1)).dof(d(1)),:),Result.U_step(Node(nid(2)).dof(d(2)),:),Result.lr_step);



