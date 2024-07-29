clear all; clc;

% =========================================================================
% Problem definition
% exact solution
% material properties
E  = 30e6;     % Young's modulus  
nu = 0.3;      % Poisson's ratio   
lambda = nu * E / ((1+nu) * (1-2 * nu));   % lame
mu = E / (2 * (1+nu)); 

D  =[lambda+2*mu,lambda,     lambda,     0, 0, 0;           
     lambda,     lambda+2*mu,lambda,     0, 0, 0;
     lambda,     lambda,     lambda+2*mu,0, 0, 0;
     0,          0,          0,          mu,0, 0;
     0,          0,          0,          0, mu,0;
     0,          0,          0,          0, 0, mu]; 

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

exact_ux   = @(x,y,z) x*(1-x)*y*(1-y)*z*(1-z) + 0.1*G(x, y, z);
exact_ux_x = @(x,y,z) (1-2*x)*y*(1-y)*z*(1-z) + 0.1*G_x(x, y, z); 
exact_ux_y = @(x,y,z) x*(1-x)*(1-2*y)*z*(1-z) + 0.1*G_y(x, y, z);
exact_ux_z = @(x,y,z) x*(1-x)*y*(1-y)*(1-2*z) + 0.1*G_z(x, y, z);
exact_ux_xx = @(x,y,z) -2*y*(1-y)*z*(1-z) + 0.1*G_xx(x, y, z); 
exact_ux_yy = @(x,y,z) -2*x*(1-x)*z*(1-z) + 0.1*G_yy(x, y, z);
exact_ux_zz = @(x,y,z) -2*x*(1-x)*y*(1-y) + 0.1*G_zz(x, y, z);
exact_ux_yz = @(x,y,z) x*(1-x)*(1-2*y)*(1-2*z) + 0.1*G_yz(x, y, z); 
exact_ux_zx = @(x,y,z) (1-2*x)*y*(1-y)*(1-2*z) + 0.1*G_zx(x, y, z); 
exact_ux_xy = @(x,y,z) (1-2*x)*(1-2*y)*z*(1-z) + 0.1*G_xy(x, y, z); 

exact_uy   = @(x,y,z) x*(1-x)*y*(1-y)*z*(1-z) + 0.1*G(x, y, z);
exact_uy_x = @(x,y,z) (1-2*x)*y*(1-y)*z*(1-z) + 0.1*G_x(x, y, z); 
exact_uy_y = @(x,y,z) x*(1-x)*(1-2*y)*z*(1-z) + 0.1*G_y(x, y, z);
exact_uy_z = @(x,y,z) x*(1-x)*y*(1-y)*(1-2*z) + 0.1*G_z(x, y, z);
exact_uy_xx = @(x,y,z) -2*y*(1-y)*z*(1-z) + 0.1*G_xx(x, y, z); 
exact_uy_yy = @(x,y,z) -2*x*(1-x)*z*(1-z) + 0.1*G_yy(x, y, z);
exact_uy_zz = @(x,y,z) -2*x*(1-x)*y*(1-y) + 0.1*G_zz(x, y, z);
exact_uy_yz = @(x,y,z) x*(1-x)*(1-2*y)*(1-2*z) + 0.1*G_yz(x, y, z); 
exact_uy_zx = @(x,y,z) (1-2*x)*y*(1-y)*(1-2*z) + 0.1*G_zx(x, y, z); 
exact_uy_xy = @(x,y,z) (1-2*x)*(1-2*y)*z*(1-z) + 0.1*G_xy(x, y, z); 

exact_uz   = @(x,y,z) x*(1-x)*y*(1-y)*z*(1-z) + 0.1*G(x, y, z);
exact_uz_x = @(x,y,z) (1-2*x)*y*(1-y)*z*(1-z) + 0.1*G_x(x, y, z); 
exact_uz_y = @(x,y,z) x*(1-x)*(1-2*y)*z*(1-z) + 0.1*G_y(x, y, z);
exact_uz_z = @(x,y,z) x*(1-x)*y*(1-y)*(1-2*z) + 0.1*G_z(x, y, z);
exact_uz_xx = @(x,y,z) -2*y*(1-y)*z*(1-z) + 0.1*G_xx(x, y, z); 
exact_uz_yy = @(x,y,z) -2*x*(1-x)*z*(1-z) + 0.1*G_yy(x, y, z);
exact_uz_zz = @(x,y,z) -2*x*(1-x)*y*(1-y) + 0.1*G_zz(x, y, z);
exact_uz_yz = @(x,y,z) x*(1-x)*(1-2*y)*(1-2*z) + 0.1*G_yz(x, y, z); 
exact_uz_zx = @(x,y,z) (1-2*x)*y*(1-y)*(1-2*z) + 0.1*G_zx(x, y, z); 
exact_uz_xy = @(x,y,z) (1-2*x)*(1-2*y)*z*(1-z) + 0.1*G_xy(x, y, z); 



% strain
exact_strainxx = @(x,y,z) exact_ux_x(x,y,z);
exact_strainyy = @(x,y,z) exact_uy_y(x,y,z);
exact_strainzz = @(x,y,z) exact_uz_z(x,y,z);
exact_strainyz = @(x,y,z) 0.5 * (exact_uy_z(x,y,z) + exact_uz_y(x,y,z));
exact_strainzx = @(x,y,z) 0.5 * (exact_uz_x(x,y,z) + exact_ux_z(x,y,z));
exact_strainxy = @(x,y,z) 0.5 * (exact_ux_y(x,y,z) + exact_uy_x(x,y,z));

% stress
exact_stressxx = @(x,y,z) (lambda + 2*mu)*exact_strainxx(x,y,z) + lambda * exact_strainyy(x,y,z)+ lambda * exact_strainzz(x,y,z);
exact_stressyy = @(x,y,z) (lambda + 2*mu)*exact_strainyy(x,y,z) + lambda * exact_strainxx(x,y,z)+ lambda * exact_strainzz(x,y,z);
exact_stresszz = @(x,y,z) (lambda + 2*mu)*exact_strainzz(x,y,z) + lambda * exact_strainxx(x,y,z)+ lambda * exact_strainyy(x,y,z);
exact_stressyz = @(x,y,z) mu*(exact_uy_z(x,y,z) + exact_uz_y(x,y,z));
exact_stresszx = @(x,y,z) mu*(exact_uz_x(x,y,z) + exact_ux_z(x,y,z));
exact_stressxy = @(x,y,z) mu*(exact_ux_y(x,y,z) + exact_uy_x(x,y,z));

exact_stressxx_x = @(x,y,z) (lambda + 2*mu)*exact_ux_xx(x,y,z) + lambda * exact_uy_xy(x,y,z) + lambda * exact_uz_zx(x,y,z) ;
exact_stressxy_y = @(x,y,z) mu*(exact_ux_yy(x,y,z) + exact_uy_xy(x,y,z)) ;
exact_stressxz_z = @(x,y,z) mu*(exact_uz_zx(x,y,z) + exact_ux_zz(x,y,z)) ;

exact_stressyy_y = @(x,y,z) (lambda + 2*mu)*exact_uy_yy(x,y,z) + lambda * exact_ux_xy(x,y,z) + lambda * exact_uz_yz(x,y,z);
exact_stressyx_x = @(x,y,z) mu*(exact_ux_xy(x,y,z) + exact_uy_xx(x,y,z)) ;
exact_stressyz_z = @(x,y,z) mu*(exact_uy_zz(x,y,z) + exact_uz_yz(x,y,z)) ;

exact_stresszz_z = @(x,y,z) (lambda + 2*mu)*exact_uz_zz(x,y,z) + lambda * exact_ux_zx(x,y,z)+ lambda * exact_uy_yz(x,y,z);
exact_stresszx_x = @(x,y,z) mu*(exact_uz_xx(x,y,z) + exact_ux_zx(x,y,z)) ;
exact_stresszy_y = @(x,y,z) mu*(exact_uy_yz(x,y,z) + exact_uz_yy(x,y,z)) ;


% force
f_x = @(x,y,z) -exact_stressxx_x(x,y,z) - exact_stressxy_y(x,y,z) - exact_stressxz_z(x,y,z);
f_y = @(x,y,z) -exact_stressyx_x(x,y,z) - exact_stressyy_y(x,y,z) - exact_stressyz_z(x,y,z);
f_z = @(x,y,z) -exact_stresszx_x(x,y,z) - exact_stresszy_y(x,y,z) - exact_stresszz_z(x,y,z);

f = {f_x, f_y, f_z};


% Dirichlet BC
g_x = @(x,y,z) exact_ux(x,y,z);
g_y = @(x,y,z) exact_uy(x,y,z);
g_z = @(x,y,z) exact_uz(x,y,z);

g = {g_x, g_y, g_z};


% Neumann BC
h1_x = @(x,y,z) -exact_stresszx(x,y,z);
h2_x = @(x,y,z) exact_stresszx(x,y,z);
h1_y = @(x,y,z) -exact_stressyz(x,y,z);
h2_y = @(x,y,z) exact_stressyz(x,y,z);
h1_z = @(x,y,z) -exact_stresszz(x,y,z);
h2_z = @(x,y,z) exact_stresszz(x,y,z);

h1 = {h1_x, h1_y, h1_z};
h2 = {h2_x, h2_y, h2_z};


% =========================================================================

% =========================================================================
% Generate the mesh
% FEM mesh settings
n_sd = 3;                 % space dimension

n_el_x = 10;               % number of element in x-direction
n_el_y = 10;               % number of element in y-direction
n_el_z = 10;               % number of element in z-direction
n_el   = n_el_x * n_el_y * n_el_z;   % total number of element in 3D domain

n_en = 8; % number of  element nodes
n_ed = 3; % number of element degrees of freedom (per node)

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
% n_int = 6;  % number of 6 point quadrature rule

% quadrature rule
n_int_xi  = 3;                                  % number of quadrature points in x-direction
n_int_eta = 3;                                  % number of quadrature points in y-direction
n_int_zeta = 3;                                 % number of quadrature points in z-direction
n_int     = n_int_xi * n_int_eta * n_int_zeta;  % number of total quadrature points

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
n_int_xi = 4; n_int_eta = 4; n_int_zeta = 4;
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
 
   x_ele = zeros(n_en, 1);         % x-coordinate of nodes in ee th element
   y_ele = zeros(n_en, 1);         % y-coordinate of nodes in ee th element
   z_ele = zeros(n_en, 1);         % z-coordinate of nodes in ee th element
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

     
     Bmat = [];
     for aa = 1 : n_en
         [Na_xi, Na_eta, Na_zeta] = PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
         Na_x = (Na_xi * conf11 + Na_eta * conf21 + Na_zeta * conf31) / detJ;
         Na_y = (Na_xi * conf12 + Na_eta * conf22 + Na_zeta * conf32) / detJ;
         Na_z = (Na_xi * conf13 + Na_eta * conf23 + Na_zeta * conf33) / detJ;

         %Set up B
         Ba = [Na_x, 0,    0;
               0,    Na_y, 0;
               0,    0,    Na_z;
               0,    Na_z, Na_y;
               Na_z, 0,    Na_x;
               Na_y, Na_x, 0];
         Bmat = [Bmat, Ba];

         % Set up f_ele
         for ii = 1 : n_ed
              pp = n_ed * (aa - 1) + ii;
              f_ele(pp) = f_ele(pp) + weight(ll) * detJ * f{ii}(x_l, y_l, z_l) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
         end
     end % end of aa-loop

     % Set up k_ele
     k_ele = k_ele + Bmat' * D * Bmat * detJ * weight(ll);

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
             for ii = 1 : n_ed
                 pp = n_ed * (aa - 1) + ii;
                 if ismember(IEN(aa, ee),h1_node_number)
                     f_ele(pp) = f_ele(pp) + weight2D(ll) * j_h1 * h1{ii}(x_l_2D_h1,y_l_2D_h1, 0) * PolyShape_3d(deg_xi,deg_eta,deg_zeta, aa, xi2D(ll), eta2D(ll), -1);
                 elseif ismember(IEN(aa, ee),h2_node_number)
                     f_ele(pp) = f_ele(pp) + weight2D(ll) * j_h2 * h2{ii}(x_l_2D_h2,y_l_2D_h2, 1) * PolyShape_3d(deg_xi,deg_eta,deg_zeta, aa, xi2D(ll), eta2D(ll), 1);
                 end
             end
        end % end of aa-loop
   end % end of quadrature loop
  

        

   % global assembly
   for aa = 1 : n_en
       for ii = 1 : n_ed
           pp = n_ed * (aa - 1) + ii;
           PP = LM(pp, ee);
           if PP > 0

               F(PP) = F(PP) + f_ele(pp);
               for bb = 1 : n_en
                   for jj = 1 : n_ed
                       qq = n_ed * (bb - 1) + jj;
                       QQ = LM(qq, ee);
                       if QQ > 0
                          K(PP, QQ) = K(PP, QQ) + k_ele(pp, qq);
                       else
                         % do something for non-zero g boundary condition
                          F(PP) = F(PP) - k_ele(pp, qq) * g{jj}(x_ele(bb), y_ele(bb), z_ele(bb));
                       end  % end of if QQ
                   end % end of jj-loop
               end   % end of for bb
           end   % end of if PP
       end % end of ii-loop
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
disp = zeros(n_np * n_ed, 1);
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
      for ii = 1: n_ed
          index  = (IEN(aa,ee) - 1) * n_ed + ii;
          pp = n_ed * (aa - 1) + ii;
          PP = LM(pp, ee);           % equation number 
          if PP > 0
              disp(index) = d_temp(PP);
          else
              disp(index) = g{ii}(x_ele(aa), y_ele(aa), z_ele(aa));
          end
      end
   end
end

% Extract the ux displacement component
ux_h = zeros(n_np, 1);
% Extract the uy displacement component
uy_h = zeros(n_np, 1);
% Extract the uz displacement component
uz_h = zeros(n_np, 1);

for ii = 1:n_np
    index = (ii - 1) * n_ed;
    ux_h(ii) = disp(index + 1);
    uy_h(ii) = disp(index + 2);
    uz_h(ii) = disp(index + 3);
end




% plot mesh
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'blue', 'EdgeColor', 'black','FaceAlpha', 0.7);

    grid on;
    view(3); % 设置 3D 视图
    hold on;
end
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');  
% title('Uneformed-Wireframe Cube');
% plot displacement
% figure;
xnew_coor = zeros(n_np, 1);
ynew_coor = zeros(n_np, 1);
znew_coor = zeros(n_np, 1);
for ii = 1 :  n_np
    xnew_coor(ii) = x_coor(ii) + ux_h(ii);
    ynew_coor(ii) = y_coor(ii) + uy_h(ii);
    znew_coor(ii) = z_coor(ii) + uz_h(ii);
end
for ee = 1 : n_el
     % 顶点坐标
    vertices = [xnew_coor(IEN(1,ee)), ynew_coor(IEN(1,ee)), znew_coor(IEN(1,ee));  % V1
                xnew_coor(IEN(2,ee)), ynew_coor(IEN(2,ee)), znew_coor(IEN(2,ee));  % V2
                xnew_coor(IEN(3,ee)), ynew_coor(IEN(3,ee)), znew_coor(IEN(3,ee));  % V3
                xnew_coor(IEN(4,ee)), ynew_coor(IEN(4,ee)), znew_coor(IEN(4,ee));  % V4
                xnew_coor(IEN(5,ee)), ynew_coor(IEN(5,ee)), znew_coor(IEN(5,ee));  % V5
                xnew_coor(IEN(6,ee)), ynew_coor(IEN(6,ee)), znew_coor(IEN(6,ee));  % V6
                xnew_coor(IEN(7,ee)), ynew_coor(IEN(7,ee)), znew_coor(IEN(7,ee));  % V7
                xnew_coor(IEN(8,ee)), ynew_coor(IEN(8,ee)), znew_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'red', 'EdgeColor', 'black','FaceAlpha', 0.7);
    
    grid on;
    view(3); % 设置 3D 视图
    hold on;
end
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');  
title('Undeformed v.s. Deformed-Wireframe Cube');



% plot the solution ux
figure;
[X_h, Y_h, Z_h] = meshgrid( 0:hh_x:1, 0:hh_y:1, 0:hh_z:1 );
UX_h= reshape(ux_h, n_np_x, n_np_y, n_np_z);
slice(X_h, Y_h, Z_h, UX_h, [], [], 0:hh_z:1);
% shading interp;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('FEM-3D X-displacement Field Slice Plot');

figure;
[X, Y, Z] = meshgrid( 0:0.01:1, 0:0.01:1, 0:0.01:1 );
for ii = 1:101
    for jj = 1:101
        for kk = 1:101
            UX_exact(ii,jj,kk) = exact_ux(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk));
        end
    end
end
slice(X,Y,Z, UX_exact, [], [], 0:hh_z:1);
% shading interp;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Exact-3D X-displacement Field Slice Plot');

% plot the solution uy
figure;
[X_h, Y_h, Z_h] = meshgrid( 0:hh_x:1, 0:hh_y:1, 0:hh_z:1 );
UY_h= reshape(uy_h, n_np_x, n_np_y, n_np_z);
slice(X_h, Y_h, Z_h, UY_h, [], [], 0:hh_z:1);
% shading interp;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('FEM-3D Y-displacement Field Slice Plot');

figure;
[X, Y, Z] = meshgrid( 0:0.01:1, 0:0.01:1, 0:0.01:1 );
for ii = 1:101
    for jj = 1:101
        for kk = 1:101
            UY_exact(ii,jj,kk) = exact_uy(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk));
        end
    end
end
slice(X,Y,Z, UY_exact, [], [], 0:hh_z:1);
% shading interp;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Exact-3D Y-displacement Field Slice Plot');

% plot the solution uz
figure;
[X_h, Y_h, Z_h] = meshgrid( 0:hh_x:1, 0:hh_y:1, 0:hh_z:1 );
UZ_h= reshape(uz_h, n_np_x, n_np_y, n_np_z);
slice(X_h, Y_h, Z_h, UZ_h, [], [], 0:hh_z:1);
% shading interp;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('FEM-3D Z-displacement Field Slice Plot');

figure;
[X, Y, Z] = meshgrid( 0:0.01:1, 0:0.01:1, 0:0.01:1 );
for ii = 1:101
    for jj = 1:101
        for kk = 1:101
            UZ_exact(ii,jj,kk) = exact_uz(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk));
        end
    end
end
slice(X,Y,Z, UZ_exact, [], [], 0:hh_z:1);
% shading interp;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Exact-3D Z-displacement Field Slice Plot');



% postprocess the solution by calculating the error measured in L2 norm
nqp = 10; % we need more points 
[xi, eta,zeta, weight] = Gauss3D(nqp, nqp, nqp);

errorL2_x = 0.0; errorL2_y = 0.0; errorL2_z = 0.0; bottomL2_x = 0.0; bottomL2_y = 0.0; bottomL2_z = 0.0;
errorH1_x = 0.0; errorH1_y = 0.0; errorH1_z = 0.0; bottomH1_x = 0.0; bottomH1_y = 0.0; bottomH1_z = 0.0;
for ee = 1 : n_el

    x_ele = x_coor( IEN(1:n_en, ee) );
    y_ele = y_coor( IEN(1:n_en, ee) );
    z_ele = z_coor( IEN(1:n_en, ee) );
    ux_ele = ux_h(  IEN(1:n_en, ee) );
    uy_ele = uy_h(  IEN(1:n_en, ee) );
    uz_ele = uz_h(  IEN(1:n_en, ee) );

    for ll = 1 : n_int
        x_l = 0.0;  y_l = 0.0;  z_l = 0.0; 
        ux_l = 0.0; uy_l = 0.0; uz_l = 0.0;
        ux_l_xi  = 0.0; ux_l_eta = 0.0; ux_l_zeta = 0.0;
        uy_l_xi  = 0.0; uy_l_eta = 0.0; uy_l_zeta = 0.0;
        uz_l_xi  = 0.0; uz_l_eta = 0.0; uz_l_zeta = 0.0;
        dx_dxi   = 0.0; dy_dxi   = 0.0; dz_dxi    = 0.0;
        dx_deta  = 0.0; dy_deta  = 0.0; dz_deta   = 0.0;
        dx_dzeta = 0.0; dy_dzeta = 0.0; dz_dzeta  = 0.0;
        for aa = 1 : n_en
            x_l  = x_l  + x_ele(aa)  * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
            y_l  = y_l  + y_ele(aa)  * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
            z_l  = z_l  + z_ele(aa)  * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
            ux_l = ux_l + ux_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
            uy_l = uy_l + uy_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
            uz_l = uz_l + uz_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));

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

            ux_l_xi   = ux_l_xi   + ux_ele(aa) * Na_xi;
            ux_l_eta  = ux_l_eta  + ux_ele(aa) * Na_eta;
            ux_l_zeta = ux_l_zeta + ux_ele(aa) * Na_zeta;

            uy_l_xi   = uy_l_xi   + uy_ele(aa) * Na_xi;
            uy_l_eta  = uy_l_eta  + uy_ele(aa) * Na_eta;
            uy_l_zeta = uy_l_zeta + uy_ele(aa) * Na_zeta;

            uz_l_xi   = uz_l_xi   + uz_ele(aa) * Na_xi;
            uz_l_eta  = uz_l_eta  + uz_ele(aa) * Na_eta;
            uz_l_zeta = uz_l_zeta + uz_ele(aa) * Na_zeta;
            
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

        ux_l_x = (ux_l_xi * conf11 + ux_l_eta * conf21 + ux_l_zeta * conf31) / detJ;
        ux_l_y = (ux_l_xi * conf12 + ux_l_eta * conf22 + ux_l_zeta * conf32) / detJ;
        ux_l_z = (ux_l_xi * conf13 + ux_l_eta * conf23 + ux_l_zeta * conf33) / detJ;

        uy_l_x = (uy_l_xi * conf11 + uy_l_eta * conf21 + uy_l_zeta * conf31) / detJ;
        uy_l_y = (uy_l_xi * conf12 + uy_l_eta * conf22 + uy_l_zeta * conf32) / detJ;
        uy_l_z = (uy_l_xi * conf13 + uy_l_eta * conf23 + uy_l_zeta * conf33) / detJ;

        uz_l_x = (uz_l_xi * conf11 + uz_l_eta * conf21 + uz_l_zeta * conf31) / detJ;
        uz_l_y = (uz_l_xi * conf12 + uz_l_eta * conf22 + uz_l_zeta * conf32) / detJ;
        uz_l_z = (uz_l_xi * conf13 + uz_l_eta * conf23 + uz_l_zeta * conf33) / detJ;


        errorL2_x = errorL2_x + weight(ll) * detJ * (ux_l - exact_ux(x_l, y_l, z_l))^2;
        errorL2_y = errorL2_y + weight(ll) * detJ * (uy_l - exact_uy(x_l, y_l, z_l))^2;
        errorL2_z = errorL2_z + weight(ll) * detJ * (uz_l - exact_uz(x_l, y_l, z_l))^2;
    
        errorH1_x = errorH1_x + weight(ll) * detJ * (( ux_l_x- exact_ux_x(x_l,y_l,z_l))^2 + ( ux_l_y - exact_ux_y(x_l,y_l,z_l))^2 + ( ux_l_z - exact_ux_z(x_l,y_l,z_l))^2 ); 
        errorH1_y = errorH1_y + weight(ll) * detJ * (( uy_l_x- exact_uy_x(x_l,y_l,z_l))^2 + ( uy_l_y - exact_uy_y(x_l,y_l,z_l))^2 + ( uy_l_z - exact_uy_z(x_l,y_l,z_l))^2 );
        errorH1_z = errorH1_z + weight(ll) * detJ * (( uz_l_x- exact_uz_x(x_l,y_l,z_l))^2 + ( uz_l_y - exact_uz_y(x_l,y_l,z_l))^2 + ( uz_l_z - exact_uz_z(x_l,y_l,z_l))^2 );
    
        bottomL2_x = bottomL2_x + weight(ll) * detJ * exact_ux(x_l, y_l, z_l)^2;
        bottomL2_y = bottomL2_y + weight(ll) * detJ * exact_uy(x_l, y_l, z_l)^2;
        bottomL2_z = bottomL2_z + weight(ll) * detJ * exact_uz(x_l, y_l, z_l)^2;
    
        bottomH1_x = bottomH1_x + weight(ll) * detJ * (exact_ux_x(x_l,y_l,z_l)^2 + exact_ux_y(x_l,y_l,z_l)^2 + exact_ux_z(x_l,y_l,z_l)^2);
        bottomH1_y = bottomH1_y + weight(ll) * detJ * (exact_uy_x(x_l,y_l,z_l)^2 + exact_uy_y(x_l,y_l,z_l)^2 + exact_uy_z(x_l,y_l,z_l)^2);
        bottomH1_z = bottomH1_z + weight(ll) * detJ * (exact_uz_x(x_l,y_l,z_l)^2 + exact_uz_y(x_l,y_l,z_l)^2 + exact_uz_z(x_l,y_l,z_l)^2);
    end
end


errorL2_x = sqrt(errorL2_x) / sqrt(bottomL2_x);
errorL2_y = sqrt(errorL2_y) / sqrt(bottomL2_y);
errorL2_z = sqrt(errorL2_z) / sqrt(bottomL2_z);
errorH1_x = sqrt(errorH1_x) / sqrt(bottomH1_x);
errorH1_y = sqrt(errorH1_y) / sqrt(bottomH1_y);
errorH1_z = sqrt(errorH1_z) / sqrt(bottomH1_z);

errorL2 = errorL2_x + errorL2_y + errorL2_z;
errorH1 = errorH1_x + errorH1_y + errorH1_z;




% get stress & stain at gauss points
n_int_xi = 4; n_int_eta = 4; n_int_zeta = 4;
n_int = n_int_xi*n_int_eta*n_int_zeta;
[xi, eta, zeta, weight] = Gauss3D(n_int_xi, n_int_eta, n_int_zeta);
for ee = 1 : n_el
    fprintf(1,'Element  %d \n',ee);
    fprintf(1,'-------------\n');

    x_ele = x_coor( IEN(1:n_en, ee) );
    y_ele = y_coor( IEN(1:n_en, ee) );
    z_ele = z_coor( IEN(1:n_en, ee) );

    disp_ele = zeros(n_ee, 1);
    strain_ele = []; stress_ele = [];
    X_ll = []; Y_ll = [];Z_ll = [];

    for ll = 1 : n_int

       x_l = 0.0;    y_l = 0.0;     z_l=0.0;                % coordinate in terms of xi(ll)
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

          % displacement at each degree of freedom of ee-th element
          for ii = 1 : n_ed
              index  = (IEN(aa,ee) - 1) * n_ed + ii;
              pp = n_ed * (aa - 1) + ii;
              disp_ele(pp) = disp(index);
          end


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

     
       Bmat = [];
       for aa = 1 : n_en
           [Na_xi, Na_eta, Na_zeta] = PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, aa, xi(ll), eta(ll), zeta(ll));
           Na_x = (Na_xi * conf11 + Na_eta * conf21 + Na_zeta * conf31) / detJ;
           Na_y = (Na_xi * conf12 + Na_eta * conf22 + Na_zeta * conf32) / detJ;
           Na_z = (Na_xi * conf13 + Na_eta * conf23 + Na_zeta * conf33) / detJ;

           %Set up B
           Ba = [Na_x, 0,    0;
                 0,    Na_y, 0;
                 0,    0,    Na_z;
                 0,    Na_z, Na_y;
                 Na_z, 0,    Na_x;
                 Na_y, Na_x, 0];
           Bmat = [Bmat, Ba];

       end % end of aa-loop

       X_ll = [X_ll, x_l];
       Y_ll = [Y_ll, y_l];
       Z_ll = [Z_ll, z_l];

       strain_ll = Bmat * disp_ele;           % strains at ll-th point
       strain_ele = [strain_ele, strain_ll];    

       stress_ll = D * strain_ll;            % stresses at ll-th point
       stress_ele = [stress_ele, stress_ll];
  
   end % end of quadrature loop


   % strains at gauss points
   strainxx = strain_ele(1, :);
   strainyy = strain_ele(2, :);
   strainzz = strain_ele(3, :);
   strainyz = strain_ele(4, :);
   strainzx = strain_ele(5, :);
   strainxy = strain_ele(6, :);
    
   % stress at gauss points
   stressxx = stress_ele(1, :);
   stressyy = stress_ele(2, :);
   stresszz = stress_ele(3, :);
   stressyz = stress_ele(4, :);
   stresszx = stress_ele(5, :);
   stressxy = stress_ele(6, :);
    
   strain_gauss = [X_ll', Y_ll', Z_ll', strainxx', strainyy', strainzz', strainyz', strainzx',strainxy'];
   fprintf('%-7s %7s %7s %7s %7s %7s %7s %7s %7s\n', 'x-coord', 'y-coord', 'z-coord', 'e_xx', 'e_yy', 'e_zz', 'e_yz', 'e_zx', 'e_xy');
   fprintf(1,'%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n',strain_gauss'); 
   
    
   stress_gauss = [X_ll', Y_ll', Z_ll', stressxx', stressyy', stresszz', stressyz', stresszx',stressxy'];
   fprintf('%-7s %7s %7s %7s %7s %7s %7s %7s %7s\n', 'x-coord', 'y-coord', 'z-coord', 's_xx', 's_yy', 's_zz', 's_yz', 's_zx', 's_xy');
   fprintf(1,'%.5f\t %.5f\t %.5f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\n',stress_gauss'); 
   
end

% get stress & strain at nodes
strain_nodes = zeros(n_np, 6);
stress_nodes = zeros(n_np, 6);
counter = zeros(n_np);

for ee = 1 : n_el
    fprintf(1,'Element  %d \n',ee);
    fprintf(1,'-------------\n');

    x_ele = x_coor( IEN(1:n_en, ee) );
    y_ele = y_coor( IEN(1:n_en, ee) );
    z_ele = z_coor( IEN(1:n_en, ee) );

    disp_ele = zeros(n_ee, 1);
    strain_ele = []; stress_ele = [];

    xi_val   = [-1,  1,  1, -1, -1, 1, 1, -1];    % xi value at nodes
    eta_val  = [-1, -1,  1,  1, -1,-1, 1,  1];    % eta value at nodes
    zeta_val = [-1, -1, -1, -1,  1, 1, 1,  1];    % zeta value at nodes

    for ll = 1 : n_en
       dx_dxi = 0.0; dx_deta = 0.0; dx_dzeta = 0.0;
       dy_dxi = 0.0; dy_deta = 0.0; dy_dzeta = 0.0;
       dz_dxi = 0.0; dz_deta = 0.0; dz_dzeta = 0.0;

       for aa = 1 : n_en
          x_l = x_l + x_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi_val(ll), eta_val(ll), zeta_val(ll));
          y_l = y_l + y_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi_val(ll), eta_val(ll), zeta_val(ll));
          z_l = z_l + z_ele(aa) * PolyShape_3d(deg_xi, deg_eta, deg_zeta, aa, xi_val(ll), eta_val(ll), zeta_val(ll));

          [Na_xi, Na_eta, Na_zeta] = PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, aa, xi_val(ll), eta_val(ll), zeta_val(ll));

          dx_dxi   = dx_dxi   + x_ele(aa) * Na_xi;
          dx_deta  = dx_deta  + x_ele(aa) * Na_eta;
          dx_dzeta = dx_dzeta + x_ele(aa) * Na_zeta;

          dy_dxi   = dy_dxi   + y_ele(aa) * Na_xi;
          dy_deta  = dy_deta  + y_ele(aa) * Na_eta;
          dy_dzeta = dy_dzeta + y_ele(aa) * Na_zeta;

          dz_dxi   = dz_dxi   + z_ele(aa) * Na_xi;
          dz_deta  = dz_deta  + z_ele(aa) * Na_eta;
          dz_dzeta = dz_dzeta + z_ele(aa) * Na_zeta;

          % displacement at each degree of freedom of ee-th element
          for ii = 1 : n_ed
              index  = (IEN(aa,ee) - 1) * n_ed + ii;
              pp = n_ed * (aa - 1) + ii;
              disp_ele(pp) = disp(index);
          end


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

     
       Bmat = [];
       for aa = 1 : n_en
           [Na_xi, Na_eta, Na_zeta] = PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, aa, xi_val(ll), eta_val(ll), zeta_val(ll));
           Na_x = (Na_xi * conf11 + Na_eta * conf21 + Na_zeta * conf31) / detJ;
           Na_y = (Na_xi * conf12 + Na_eta * conf22 + Na_zeta * conf32) / detJ;
           Na_z = (Na_xi * conf13 + Na_eta * conf23 + Na_zeta * conf33) / detJ;

           %Set up B
           Ba = [Na_x, 0,    0;
                 0,    Na_y, 0;
                 0,    0,    Na_z;
                 0,    Na_z, Na_y;
                 Na_z, 0,    Na_x;
                 Na_y, Na_x, 0];
           Bmat = [Bmat, Ba];

       end % end of aa-loop

       strain_ll = Bmat * disp_ele;           % strains at ll-th point
       strain_ele = [strain_ele, strain_ll];    

       stress_ll = D * strain_ll;            % stresses at ll-th point
       stress_ele = [stress_ele, stress_ll];
  
   end % end of quadrature loop


   % strains at gauss points
   strainxx = strain_ele(1, :);
   strainyy = strain_ele(2, :);
   strainzz = strain_ele(3, :);
   strainyz = strain_ele(4, :);
   strainzx = strain_ele(5, :);
   strainxy = strain_ele(6, :);
    
   % stress at gauss points
   stressxx = stress_ele(1, :);
   stressyy = stress_ele(2, :);
   stresszz = stress_ele(3, :);
   stressyz = stress_ele(4, :);
   stresszx = stress_ele(5, :);
   stressxy = stress_ele(6, :);
    
   strain_node = [x_ele, y_ele, z_ele, strainxx', strainyy', strainzz', strainyz', strainzx',strainxy'];
   fprintf('%-7s %7s %7s %7s %7s %7s %7s %7s %7s\n', 'x-coord', 'y-coord', 'z-coord', 'e_xx', 'e_yy', 'e_zz', 'e_yz', 'e_zx', 'e_xy');
   fprintf(1,'%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n',strain_node'); 
   
    
   stress_node = [x_ele, y_ele, z_ele, stressxx', stressyy', stresszz', stressyz', stresszx',stressxy'];
   fprintf('%-7s %7s %7s %7s %7s %7s %7s %7s %7s\n', 'x-coord', 'y-coord', 'z-coord', 's_xx', 's_yy', 's_zz', 's_yz', 's_zx', 's_xy');
   fprintf(1,'%.5f\t %.5f\t %.5f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\n',stress_node'); 

   % essemble the element strain & stress to global strain & stress
   strain_nodes(IEN(:,ee), :) = strain_node(:, 4:9);                
   stress_nodes(IEN(:,ee), :) = stress_node(:, 4:9);

   counter( IEN(:,ee) ) = counter( IEN(:, ee) ) + ones(n_en,1);    % count the time a stress is added to a node 
   
end

% plot stress and strain contour
% exx------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    exx = strain_nodes(IEN(:,ee),1)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', exx, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\epsilon_x_x - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% eyy------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    eyy = strain_nodes(IEN(:,ee),2)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', eyy, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\epsilon_y_y - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% ezz------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    ezz = strain_nodes(IEN(:,ee),3)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', ezz, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\epsilon_z_z - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% eyz------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    eyz = strain_nodes(IEN(:,ee),4)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', eyz, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\epsilon_y_z - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% ezx------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    ezx = strain_nodes(IEN(:,ee),5)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', ezx, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\epsilon_z_x - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% eyy------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    exy = strain_nodes(IEN(:,ee),6)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', exy, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\epsilon_x_y - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% sxx------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    sxx = stress_nodes(IEN(:,ee),1)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', sxx, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\sigma_x_x - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% syy------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    syy = stress_nodes(IEN(:,ee),2)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', syy, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\sigma_y_y - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% szz------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    szz = stress_nodes(IEN(:,ee),3)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', szz, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\sigma_z_z - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% syz------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    syz = stress_nodes(IEN(:,ee),4)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', syz, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\sigma_y_z - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% szx------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    szx = stress_nodes(IEN(:,ee),5)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', szx, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\sigma_z_x - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% sxy------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    sxy = stress_nodes(IEN(:,ee),6)./counter( IEN(:,ee) ); 

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', sxy, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('\sigma_x_y - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;

% Von Mises------------------------------------------------------------
figure;
for ee = 1 : n_el
    % 顶点坐标
    vertices = [x_coor(IEN(1,ee)), y_coor(IEN(1,ee)), z_coor(IEN(1,ee));  % V1
                x_coor(IEN(2,ee)), y_coor(IEN(2,ee)), z_coor(IEN(2,ee));  % V2
                x_coor(IEN(3,ee)), y_coor(IEN(3,ee)), z_coor(IEN(3,ee));  % V3
                x_coor(IEN(4,ee)), y_coor(IEN(4,ee)), z_coor(IEN(4,ee));  % V4
                x_coor(IEN(5,ee)), y_coor(IEN(5,ee)), z_coor(IEN(5,ee));  % V5
                x_coor(IEN(6,ee)), y_coor(IEN(6,ee)), z_coor(IEN(6,ee));  % V6
                x_coor(IEN(7,ee)), y_coor(IEN(7,ee)), z_coor(IEN(7,ee));  % V7
                x_coor(IEN(8,ee)), y_coor(IEN(8,ee)), z_coor(IEN(8,ee));  % V8
                ];
    % 面
    faces = [1,2,3,4;  % 底面
             5,6,7,8;  % 顶面
             1,2,6,5;  % 侧面1
             2,3,7,6;  % 侧面2
             3,4,8,7;  % 侧面3
             4,1,5,8;  % 侧面4
             ];

    sxx = stress_nodes(IEN(:,ee),1)./counter( IEN(:,ee) ); 
    syy = stress_nodes(IEN(:,ee),2)./counter( IEN(:,ee) );
    szz = stress_nodes(IEN(:,ee),3)./counter( IEN(:,ee) );
    syz = stress_nodes(IEN(:,ee),4)./counter( IEN(:,ee) );
    szx = stress_nodes(IEN(:,ee),5)./counter( IEN(:,ee) );
    sxy = stress_nodes(IEN(:,ee),6)./counter( IEN(:,ee) );

   
    mises = sqrt(0.5 * ((sxx - syy).^2 + (syy - szz).^2 + (szz - sxx).^2) + 3 * (sxy.^2 + syz.^2 + szx.^2));

    % 创建 patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', mises, 'FaceColor', 'interp', 'EdgeColor', 'black');

    hold on;
end
% 添加颜色条
colorbar;
% 设置轴标签
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Von Mises \sigma - 3D Cloud Plot with Vertex Values');
% 设置视图和比例
view(3);
axis equal;
grid on;




% EOF