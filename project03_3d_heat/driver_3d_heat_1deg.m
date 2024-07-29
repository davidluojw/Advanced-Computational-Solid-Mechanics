clear all; clc;

% =========================================================================
% Problem definition
% exact solution
% manufactured solution and source term
G = @(x, y, z) sin((x + y+ z) * 2 * pi);
G_x = @(x, y, z) 2 * pi * cos((x + y + z) * 2 * pi);
G_y = @(x, y, z) 2 * pi * cos((x + y + z) * 2 * pi);
G_z = @(x, y, z) 2 * pi * cos((x + y + z) * 2 * pi);
G_xx = @(x, y, z) -4 * pi * pi * sin((x + y + z) * 2 * pi);
G_yy = @(x, y, z) -4 * pi * pi * sin((x + y + z) * 2 * pi);
G_zz = @(x, y, z) -4 * pi * pi * sin((x + y + z) * 2 * pi);
G_yz = @(x, y, z) -4 * pi * pi * sin((x + y + z) * 2 * pi);
G_zx = @(x, y, z) -4 * pi * pi * sin((x + y + z) * 2 * pi);
G_xy = @(x, y, z) -4 * pi * pi * sin((x + y + z) * 2 * pi);

exact   = @(x,y,z) x*(1-x)*y*(1-y)*z*(1-z) + 0.1*G(x, y, z);
exact_x = @(x,y,z) (1-2*x)*y*(1-y)*z*(1-z) + 0.1*G_x(x, y, z); 
exact_y = @(x,y,z) x*(1-x)*(1-2*y)*z*(1-z) + 0.1*G_y(x, y, z);
exact_z = @(x,y,z) x*(1-x)*y*(1-y)*(1-2*z) + 0.1*G_z(x, y, z);
exact_xx = @(x,y,z) -2*y*(1-y)*z*(1-z) + 0.1*G_xx(x, y, z); 
exact_yy = @(x,y,z) -2*x*(1-x)*z*(1-z) + 0.1*G_yy(x, y, z);
exact_zz = @(x,y,z) -2*x*(1-x)*y*(1-y) + 0.1*G_zz(x, y, z);
exact_yz = @(x,y,z) x*(1-x)*(1-2*y)*(1-2*z) + 0.1*G_yz(x, y, z); 
exact_zx = @(x,y,z) (1-2*x)*y*(1-y)*(1-2*z) + 0.1*G_zx(x, y, z); 
exact_xy = @(x,y,z) (1-2*x)*(1-2*y)*z*(1-z) + 0.1*G_xy(x, y, z); 

kappa = eye(3); % isotropic homogeneous heat conductivity

% force
f = @(x,y,z) -kappa(1,1) * exact_xx(x,y,z) - kappa(2,2) * exact_yy(x,y,z)- kappa(3,3) * exact_zz(x,y,z) ... 
           -2 * kappa(2,3) * exact_yz(x,y,z) - 2 * kappa(3,1) * exact_zx(x,y,z) - 2 * kappa(1,2) * exact_xy(x,y,z) ;

% Dirichlet BC
g = @(x,y,z) 0.1*G(x, y, z);

% Neumann BC
h1 = @(x,y,z) -kappa(3,1)*exact_x(x,y,z) - kappa(3,2) * exact_y(x,y,z) - kappa(3,3)*exact_z(x,y,z);
h2 = @(x,y,z) kappa(3,1)*exact_x(x,y,z) + kappa(3,2) * exact_y(x,y,z) + kappa(3,3)*exact_z(x,y,z);


% =========================================================================

% =========================================================================
% Generate the mesh
% FEM mesh settings
n_sd = 3;                 % space dimension

n_el_x = 32;               % number of element in x-direction
n_el_y = 32;               % number of element in y-direction
n_el_z = 32;               % number of element in z-direction
n_el   = n_el_x * n_el_y * n_el_z;   % total number of element in 3D domain

n_en = 8; % number of  element nodes
n_ed = 1; % number of element degrees of freedom (per node)

deg_xi = 1;   % degree of lagrange polynomial in xi-direction
deg_eta = 1;  % degree of lagrange polynomial in eta-direction
deg_zeta = 1;  % degree of lagrange polynomial in zeta-direction

n_np_x = n_el_x * deg_xi  + 1;      % number of node points in x-direction
n_np_y = n_el_y * deg_eta + 1;      % number of node points in y-direction
n_np_z = n_el_z * deg_zeta + 1;     % number of node points in z-direction
n_np   = n_np_x * n_np_y * n_np_z;  % total number of node points in 3D domain

n_eq = (n_np - 2 * n_np_x * n_np_z - 2 * n_np_y * n_np_z + 4*n_np_z) * n_ed;   % number of equations

n_ee = n_ed * n_en;     % number of element equations

% quadrature rule
n_int = 6;  % number of 6 point quadrature rule

% =========================================================================
% generate the coordinates of the nodal points
x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);
z_coor = zeros(n_np, 1);

hh_x = 1 / (n_el_x * deg_xi);         % mesh size in the x-direction
hh_y = 1 / (n_el_y * deg_eta);        % mesh size in the y-direction
hh_z = 1 / (n_el_z * deg_zeta);        % mesh size in the z-direction

for nz = 1 : n_np_z
    for ny = 1 : n_np_y
        for nx = 1 : n_np_x
            index = (nz-1)*n_np_x*n_np_y + (ny-1)*n_np_x + nx;   % nodal index for the (nx, ny, nz) node
            x_coor(index) = (nx-1) * hh_x;
            y_coor(index) = (ny-1) * hh_y;
            z_coor(index) = (nz-1) * hh_z;
        end
    end
end

% ID and LM arrays are generated based on the BC info
ID = zeros(n_ed, n_np);
counter = 1;
for nz = 1 : n_np_z
    for ny = 1 : n_np_y
        for nx = 1:n_np_x
            AA = (nz - 1)*n_np_x * n_np_y + (ny - 1) * n_np_x + nx;
    
            % Modify ID according to the Dirichlet BC info
            if nx ~= 1 &&  nx ~= n_np_x && ny ~= 1 && ny ~= n_np_y
                for ii = 1 : n_ed
                    ID(ii, AA) = counter;
                    counter = counter + 1;
                end
            end
        end

    end
end


% trilinear Hexahedral element
% setup the IEN array for element with local node numbering as
%    a=8---------a=7
%   / |         /|
%  /  |        / |
% /   |       /  |
% a=5-|------a=6 |
% |   a=4----|---a=3
% |  /       |  /
% | /        | / 
% |/         |/
% a=1--------a=2
IEN = zeros(n_en, n_el);
for ez = 1 : n_el_z
    for ey = 1 : n_el_y
        for ex = 1:n_el_x
            ee = (ez-1)*n_el_x*n_el_y + (ey-1)*n_el_x + ex;
            IEN(1,ee) = (deg_zeta * ez - 1) * n_np_x * n_np_y + (deg_eta * ey - 1)* n_np_x + deg_xi * (ex - 1) + 1;
            IEN(2,ee) = (deg_zeta * ez - 1) * n_np_x * n_np_y + (deg_eta * ey - 1)* n_np_x + deg_xi * (ex - 1) + 2;
            IEN(3,ee) = (deg_zeta * ez - 1) * n_np_x * n_np_y +  deg_eta * ey     * n_np_x + deg_xi * (ex - 1) + 2;
            IEN(4,ee) = (deg_zeta * ez - 1) * n_np_x * n_np_y +  deg_eta * ey     * n_np_x + deg_xi * (ex - 1) + 1;
            IEN(5,ee) =  deg_zeta * ez      * n_np_x * n_np_y + (deg_eta * ey - 1)* n_np_x + deg_xi * (ex - 1) + 1;
            IEN(6,ee) =  deg_zeta * ez      * n_np_x * n_np_y + (deg_eta * ey - 1)* n_np_x + deg_xi * (ex - 1) + 2;
            IEN(7,ee) =  deg_zeta * ez      * n_np_x * n_np_y +  deg_eta * ey     * n_np_x + deg_xi * (ex - 1) + 2;
            IEN(8,ee) =  deg_zeta * ez      * n_np_x * n_np_y +  deg_eta * ey     * n_np_x + deg_xi * (ex - 1) + 1;
        end
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

% Neumann B.C. node number
h1_node_number = [];
h2_node_number = [];
for nz = 1 : n_np_z
    for ny = 1 : n_np_y
        for nx = 1:n_np_x
            AA = (nz - 1)*n_np_x*n_np_y + (ny - 1) * n_np_x + nx;
    
            % elements corresponding to the Neumann BC
            if nz == 1
                h1_node_number(end+1) = AA;
            elseif nz == n_np_z
                h2_node_number(end+1) = AA;
            end
        end
    end
end

% =========================================================================
% generate the quadrature rule
% [xi, eta, zeta, weight] = SixPts3D();
n_int_xi = 3; n_int_eta = 3; n_int_zeta = 3;
n_int = n_int_xi*n_int_eta*n_int_zeta;
[xi, eta, zeta, weight] = Gauss3D(n_int_xi, n_int_eta, n_int_zeta);
[xi2D, eta2D, weight2D] = Gauss2D(n_int_xi, n_int_eta);
% =========================================================================

% Start the assembly procedure
K = spalloc(n_eq, n_eq, 81*n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
   k_ele = zeros(n_ee, n_ee);      % dimension of element stiffness is n_ee x n_ee
   f_ele = zeros(n_ee, 1);         % element force nodes
 
   x_ele = zeros(n_en, 1);         % coordinate of nodes in ee th element
   y_ele = x_ele;
   z_ele = x_ele;
   for aa = 1 : n_en
     x_ele(aa) = x_coor( IEN(aa,ee) );  % from element node number to global node number 
     y_ele(aa) = y_coor( IEN(aa,ee) );  % A = IEN(a,e)
     z_ele(aa) = z_coor( IEN(aa,ee) );
   end

   % loop over quadrature points   
   for ll = 1 : n_int
     x_l = 0.0; y_l = 0.0;z_l=0.0;                % coordinate in terms of xi(ll)
     dx_dxi = 0.0; dx_deta = 0.0; dx_dzeta = 0.0;
     dy_dxi = 0.0; dy_deta = 0.0; dy_dzeta = 0.0;
     dz_dxi = 0.0; dz_deta = 0.0; dz_dzeta = 0.0;
     for aa = 1 : n_en
        x_l = x_l + x_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
        y_l = y_l + y_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
        z_l = z_l + z_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
        [Na_xi, Na_eta, Na_zeta] = PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
        dx_dxi   = dx_dxi   + x_ele(aa) * Na_xi;
        dx_deta  = dx_deta  + x_ele(aa) * Na_eta;
        dx_dzeta = dx_dzeta + x_ele(aa) * Na_zeta;
        dy_dxi   = dy_dxi   + y_ele(aa) * Na_xi;
        dy_deta  = dy_deta  + y_ele(aa) * Na_eta;
        dy_dzeta = dy_dzeta + y_ele(aa) * Na_zeta;
        dz_dxi   = dz_dxi   + z_ele(aa) * Na_xi;
        dz_deta  = dz_deta  + z_ele(aa) * Na_eta;
        dz_dzeta = dz_dzeta + z_ele(aa) * Na_zeta;
     end

     detJ = dx_dxi  * (dy_deta * dz_dzeta - dy_dzeta * dz_deta) ...
          - dx_deta * (dy_dxi  * dz_dzeta - dy_dzeta * dz_dxi ) ...
          + dx_dzeta* (dy_dxi  * dz_deta  - dy_deta  * dz_dxi );

     conf11 = dy_deta  * dz_dzeta  - dy_dzeta * dz_deta;
     conf12 = dy_dzeta * dz_dxi    - dy_dxi   * dz_dzeta;
     conf13 = dy_dxi   * dz_deta   - dy_deta  * dz_dxi;
     conf21 = dx_dzeta * dz_deta   - dx_deta  * dz_dzeta;
     conf22 = dx_dxi   * dz_dzeta  - dx_dzeta * dz_dxi;
     conf23 = dx_deta  * dz_dxi    - dx_dxi   * dz_deta;
     conf31 = dx_deta  * dy_dzeta  - dx_dzeta * dy_deta;
     conf32 = dx_dzeta * dy_dxi    - dx_dxi   * dy_dzeta;
     conf33 = dx_dxi   * dy_deta   - dx_deta  * dy_dxi;

     
     for aa = 1 : n_en
       [Na_xi, Na_eta, Na_zeta] = PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
       Na_x = (Na_xi * conf11 + Na_eta * conf21 + Na_zeta * conf31) / detJ;
       Na_y = (Na_xi * conf12 + Na_eta * conf22 + Na_zeta * conf32) / detJ;
       Na_z = (Na_xi * conf13 + Na_eta * conf23 + Na_zeta * conf33) / detJ;

       f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l, z_l) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
       for bb = 1 : n_en
         [Nb_xi, Nb_eta, Nb_zeta] = PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, bb, xi(ll), eta(ll), zeta(ll));
         Nb_x = (Nb_xi * conf11 + Nb_eta * conf21 + Nb_zeta * conf31) / detJ;
         Nb_y = (Nb_xi * conf12 + Nb_eta * conf22 + Nb_zeta * conf32) / detJ;
         Nb_z = (Nb_xi * conf13 + Nb_eta * conf23 + Nb_zeta * conf33) / detJ;

         k_ele(aa,bb) = k_ele(aa,bb) + weight(ll) * detJ * [Na_x, Na_y, Na_z] * kappa * [Nb_x; Nb_y;Nb_z];
       end % end of bb-loop
     end % end of aa-loop
   end % end of quadrature loop
   
   
   % loop over quadrature points for h boundary condtition 
   for ll = 1:n_int_xi * n_int_eta
        x_l_2D_h1 = 0.0; y_l_2D_h1 = 0.0;
        x_l_2D_h2 = 0.0; y_l_2D_h2 = 0.0; 
        dx_dxi_2D_h1 = 0.0; dx_deta_2D_h1 = 0.0; 
        dy_dxi_2D_h1 = 0.0; dy_deta_2D_h1 = 0.0;
        dx_dxi_2D_h2 = 0.0; dx_deta_2D_h2 = 0.0;
        dy_dxi_2D_h2 = 0.0; dy_deta_2D_h2 = 0.0;
        for aa = 1: n_en
            if ismember(IEN(aa, ee),h1_node_number)
                x_l_2D_h1 = x_l_2D_h1 + x_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi2D(ll), eta2D(ll), -1);
                y_l_2D_h1 = y_l_2D_h1 + y_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi2D(ll), eta2D(ll), -1);

                [Na_xi, Na_eta, Na_zeta] = PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, aa, xi2D(ll), eta2D(ll), -1);
                dx_dxi_2D_h1   = dx_dxi_2D_h1  + x_ele(aa) * Na_xi;
                dx_deta_2D_h1  = dx_deta_2D_h1 + x_ele(aa) * Na_eta;
                dy_dxi_2D_h1   = dy_dxi_2D_h1  + y_ele(aa) * Na_xi;
                dy_deta_2D_h1  = dy_deta_2D_h1 + y_ele(aa) * Na_eta;

            elseif ismember(IEN(aa, ee),h2_node_number)
                x_l_2D_h2 = x_l_2D_h2 + x_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi2D(ll), eta2D(ll), 1);
                y_l_2D_h2 = y_l_2D_h2 + y_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi2D(ll), eta2D(ll), 1);

                [Na_xi, Na_eta, Na_zeta] = PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, aa, xi2D(ll), eta2D(ll), 1);
                dx_dxi_2D_h2   = dx_dxi_2D_h2  + x_ele(aa) * Na_xi;
                dx_deta_2D_h2  = dx_deta_2D_h2 + x_ele(aa) * Na_eta;
                dy_dxi_2D_h2   = dy_dxi_2D_h2  + y_ele(aa) * Na_xi;
                dy_deta_2D_h2  = dy_deta_2D_h2 + y_ele(aa) * Na_eta;
            end
        end
        j_h1 = dx_dxi_2D_h1 * dy_deta_2D_h1 - dx_deta_2D_h1 * dy_dxi_2D_h1;
        j_h2 = dx_dxi_2D_h2 * dy_deta_2D_h2 - dx_deta_2D_h2 * dy_dxi_2D_h2;

        for aa = 1:n_en
           if ismember(IEN(aa, ee),h1_node_number)
               f_ele(aa) = f_ele(aa) + weight2D(ll) * j_h1 * h1(x_l_2D_h1,y_l_2D_h1, 0) * PolyShape_3d(deg_xi,deg_eta,deg_zeta, aa, xi2D(ll), eta2D(ll), -1);

           elseif ismember(IEN(aa, ee),h2_node_number)
               f_ele(aa) = f_ele(aa) + weight2D(ll) * j_h2 * h2(x_l_2D_h2,y_l_2D_h2, 1) * PolyShape_3d(deg_xi,deg_eta,deg_zeta, aa, xi2D(ll), eta2D(ll), 1);

           end
        end % end of aa-loop
   end % end of quadrature loop
  

        

   % global assembly
   for aa = 1 : n_en
     PP = LM(aa, ee);
     if PP > 0
       F(PP) = F(PP) + f_ele(aa);
       for bb = 1 : n_en
         QQ = LM(bb, ee);
         if QQ > 0
           K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
         else
           % do something for non-zero g boundary condition
           F(PP) = F(PP) - k_ele(aa, bb) * g(x_ele(bb), y_ele(bb), z_ele(bb));
         end % end of if QQ
       end   % end of for bb
       

     end   % end of if PP
   end    % end of for aa
end % end of element loop

% =========================================================================
% Now we have K and F assembled and we solve the linear system Kd = F
d_temp = K \ F;

% Solve linear system using gmres x = gmres(A,b,restart,tol)
% restart = 10000;
% tol = 10e-6;
% d_temp = gmres(K,F,restart,tol);
% =========================================================================

% Generate the full solution vector by inserting back the Dirichlet value
disp = zeros(n_np, 1);
for ee = 1: n_el
  x_ele = zeros(n_en, 1);
  y_ele = x_ele;
  z_ele = x_ele;
  for aa = 1 : n_en
     x_ele(aa) = x_coor( IEN(aa,ee) );
     y_ele(aa) = y_coor( IEN(aa,ee) );
     z_ele(aa) = z_coor( IEN(aa,ee) );
  end
  
  for aa = 1:n_en
    index = LM(aa, ee);           % equation number 
    AA = IEN(aa,ee);              % global node number 
    if index > 0
        disp(AA) = d_temp(index);
    else
        disp(AA) = g(x_ele(aa), y_ele(aa), z_ele(aa));
    end
  end
end

% plot the solution
figure;
[X_h, Y_h, Z_h] = meshgrid( 0:hh_x:1, 0:hh_y:1, 0:hh_z:1 );
temp_h= reshape(disp, n_np_x, n_np_y, n_np_z);
slice(X_h, Y_h, Z_h, temp_h, [], [], 0:hh_z:1);
% shading interp;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('FEM-3D Temperature Field Slice Plot');

figure;
[X, Y, Z] = meshgrid( 0:0.01:1, 0:0.01:1, 0:0.01:1 );
for ii = 1:101
    for jj = 1:101
        for kk = 1:101
            temp_exact(ii,jj,kk) = exact(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk));
        end
    end
end
slice(X,Y,Z, temp_exact, [], [], 0:hh_z:1);
% shading interp;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Exact-3D Temperature Field Slice Plot');




% postprocess the solution by calculating the error measured in L2 norm
nqp = 10; % we need more points 
[xi, eta,zeta, weight] = Gauss3D(nqp, nqp, nqp);

errorL2 = 0.0; bottomL2 = 0.0;
errorH1 = 0.0; bottomH1 = 0.0;
for ee = 1 : n_el

    x_ele = x_coor( IEN(1:n_en, ee) );
    y_ele = y_coor( IEN(1:n_en, ee) );
    z_ele = z_coor( IEN(1:n_en, ee) );
    u_ele = disp(   IEN(1:n_en, ee) );

    for ll = 1 : n_int
        x_l = 0.0; y_l = 0.0; z_l = 0.0; 
        u_l = 0.0;
        u_l_xi   = 0.0; u_l_eta  = 0.0; u_l_zeta = 0.0;
        dx_dxi   = 0.0; dy_dxi   = 0.0; dz_dxi   = 0.0;
        dx_deta  = 0.0; dy_deta  = 0.0; dz_deta  = 0.0;
        dx_dzeta = 0.0; dy_dzeta = 0.0; dz_dzeta = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
            y_l = y_l + y_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
            z_l = z_l + z_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
            u_l = u_l + u_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
            [Na_xi, Na_eta, Na_zeta] = PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
            dx_dxi   = dx_dxi   + x_ele(aa) * Na_xi;
            dx_deta  = dx_deta  + x_ele(aa) * Na_eta;
            dx_dzeta = dx_dzeta + x_ele(aa) * Na_zeta;
            dy_dxi   = dy_dxi   + y_ele(aa) * Na_xi;
            dy_deta  = dy_deta  + y_ele(aa) * Na_eta;
            dy_dzeta = dy_dzeta + y_ele(aa) * Na_zeta;
            dz_dxi   = dz_dxi   + z_ele(aa) * Na_xi;
            dz_deta  = dz_deta  + z_ele(aa) * Na_eta;
            dz_dzeta = dz_dzeta + z_ele(aa) * Na_zeta;
            u_l_xi   = u_l_xi   + u_ele(aa) * Na_xi;
            u_l_eta  = u_l_eta  + u_ele(aa) * Na_eta;
            u_l_zeta = u_l_zeta + u_ele(aa) * Na_zeta;
        end
        detJ = dx_dxi  * (dy_deta * dz_dzeta - dy_dzeta * dz_deta) ...
              - dx_deta * (dy_dxi  * dz_dzeta - dy_dzeta * dz_dxi ) ...
              + dx_dzeta* (dy_dxi  * dz_deta  - dy_deta  * dz_dxi );

        conf11 = dy_deta  * dz_dzeta  - dy_dzeta * dz_deta;
        conf12 = dy_dzeta * dz_dxi    - dy_dxi   * dz_dzeta;
        conf13 = dy_dxi   * dz_deta   - dy_deta  * dz_dxi;
        conf21 = dx_dzeta * dz_deta   - dx_deta  * dz_dzeta;
        conf22 = dx_dxi   * dz_dzeta  - dx_dzeta * dz_dxi;
        conf23 = dx_deta  * dz_dxi    - dx_dxi   * dz_deta;
        conf31 = dx_deta  * dy_dzeta  - dx_dzeta * dy_deta;
        conf32 = dx_dzeta * dy_dxi    - dx_dxi   * dy_dzeta;
        conf33 = dx_dxi   * dy_deta   - dx_deta  * dy_dxi;

        u_l_x = (u_l_xi * conf11 + u_l_eta * conf21 + u_l_zeta * conf31) / detJ;
        u_l_y = (u_l_xi * conf12 + u_l_eta * conf22 + u_l_zeta * conf32) / detJ;
        u_l_z = (u_l_xi * conf13 + u_l_eta * conf23 + u_l_zeta * conf33) / detJ;


        errorL2 = errorL2 + weight(ll) * detJ * ( u_l   - exact(x_l, y_l, z_l))^2;
        errorH1 = errorH1 + weight(ll) * detJ *(( u_l_x - exact_x(x_l,y_l,z_l))^2 ... 
                                              + ( u_l_y - exact_y(x_l,y_l,z_l))^2 ...
                                              + ( u_l_z - exact_z(x_l,y_l,z_l))^2);
        
        bottomL2 = bottomL2 + weight(ll) * detJ *  exact(x_l, y_l, z_l)^2;
        bottomH1 = bottomH1 + weight(ll) * detJ * (exact_x(x_l,y_l,z_l)^2 ...
                                                 + exact_y(x_l,y_l,z_l)^2 ...
                                                 + exact_z(x_l,y_l,z_l)^2);
    end
end

errorL2 = sqrt(errorL2) / sqrt(bottomL2);
errorH1 = sqrt(errorH1) / sqrt(bottomH1);

% EOF