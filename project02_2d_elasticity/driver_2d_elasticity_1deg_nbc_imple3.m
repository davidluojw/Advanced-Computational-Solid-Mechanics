clear all; clc;

% =========================================================================
% Problem definition
% exact solution
% material properties
E  = 30e6;     % Young's modulus  
nu = 0.3;      % Poisson's ratio   
lambda = nu * E / ((1+nu) * (1-2 * nu));   % lame
mu = E / (2 * (1+nu)); 
lambda_bar  = 2 * lambda * mu / (lambda + 2 * mu);
D  =[lambda_bar + 2 * mu    lambda_bar              0           
     lambda_bar             lambda_bar + 2 * mu     0   
     0                      0                       mu]; 


% geometry 
heit = 0.5;
len = 1;

% exact solution
% manufactured solution and source term
G    = @(x, y) sin((x + y) * 2 * pi);
G_x  = @(x, y) 2 * pi * cos((x + y) * 2 * pi);
G_y  = @(x, y) 2 * pi * cos((x + y) * 2 * pi);
G_xx = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);
G_yy = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);
G_xy = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);
G_yx = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);

exact_ux    = @(x,y) x*(len-x)*y*(1-y) + 0.1*G(x, y);
exact_ux_x  = @(x,y) (len-2*x)*y*(1-y) + 0.1*G_x(x, y); 
exact_ux_y  = @(x,y) x*(len-x)*(1-2*y) + 0.1*G_y(x, y);
exact_ux_xx = @(x,y) -2*y*(1-y)        + 0.1*G_xx(x, y); 
exact_ux_xy = @(x,y) (len-2*x)*(1-2*y) + 0.1*G_xy(x, y); 
exact_ux_yx = @(x,y) (len-2*x)*(1-2*y) + 0.1*G_yx(x, y);
exact_ux_yy = @(x,y) -2*x*(len-x)      + 0.1*G_yy(x, y);

exact_uy    = @(x,y) x*(len-x)*y*(1-y) + 0.1*G(x, y);
exact_uy_x  = @(x,y) (len-2*x)*y*(1-y) + 0.1*G_x(x, y);
exact_uy_y  = @(x,y) x*(len-x)*(1-2*y) + 0.1*G_y(x, y);
exact_uy_xx = @(x,y) -2*y*(1-y)        + 0.1*G_xx(x, y);
exact_uy_xy = @(x,y) (len-2*x)*(1-2*y) + 0.1*G_xy(x, y);
exact_uy_yx = @(x,y) (len-2*x)*(1-2*y) + 0.1*G_yx(x, y);
exact_uy_yy = @(x,y) -2*x*(len-x)      + 0.1*G_yy(x, y);



% strain
exact_strainxx = @(x,y) exact_ux_x(x,y);
exact_strainyy = @(x,y) exact_uy_y(x,y);
exact_strainxy = @(x,y) 0.5 * (exact_ux_y(x,y) + exact_uy_x(x,y));

% stress
exact_stressxx = @(x,y) (lambda_bar + 2*mu)*exact_strainxx(x,y) + lambda_bar * exact_strainyy(x,y);
exact_stressyy = @(x,y) (lambda_bar + 2*mu)*exact_strainyy(x,y) + lambda_bar * exact_strainxx(x,y);
exact_stressxy = @(x,y) 2*mu*exact_strainxy(x,y);

exact_stressxx_x = @(x,y) (lambda_bar + 2*mu)*exact_ux_xx(x,y) + lambda_bar * exact_uy_yx(x,y);
exact_stressyy_y = @(x,y) (lambda_bar + 2*mu)*exact_uy_yy(x,y) + lambda_bar * exact_ux_xy(x,y);
exact_stressxy_x = @(x,y) mu*(exact_ux_yx(x,y) + exact_uy_xx(x,y));
exact_stressxy_y = @(x,y) mu*(exact_ux_yy(x,y) + exact_uy_xy(x,y));


% force
f_x = @(x,y) -exact_stressxx_x(x,y) - exact_stressxy_y(x,y);
f_y = @(x,y) -exact_stressxy_x(x,y) - exact_stressyy_y(x,y);

f = {f_x, f_y};


% Dirichlet BC
g_x = @(x, y) exact_ux(x,y);
g_y = @(x,y) exact_uy(x,y);

g = {g_x, g_y};


% Neumann BC
h1_x = @(x, y) exact_stressxy(x,y);
h2_x = @(x, y) sqrt(2)/2*(exact_stressxx(x,y) - exact_stressxy(x,y));
h1_y = @(x, y) exact_stressyy(x,y);
h2_y = @(x, y) sqrt(2)/2*(exact_stressxy(x,y) - exact_stressyy(x,y));

h1 = {h1_x, h1_y};
h2 = {h2_x, h2_y};


% =========================================================================

% =========================================================================
% Generate the mesh
% FEM mesh settings
n_sd = 2;                 % space dimension

n_el_x = 8;               % number of element in x-direction
n_el_y = 4;               % number of element in y-direction
n_el   = n_el_x * n_el_y;   % total number of element in 2D domain

n_en = 4; % number of  element nodes
n_ed = 2; % number of element degrees of freedom (per node)

deg_xi = 1;   % degree of lagrange polynomial in xi-direction
deg_eta = 1;  % degree of lagrange polynomial in eta-direction

n_np_x = n_el_x * deg_xi  + 1;      % number of node points in x-direction
n_np_y = n_el_y * deg_eta + 1;      % number of node points in y-direction
n_np   = n_np_x * n_np_y; % total number of node points in 2D domain

n_eq = (n_np - 2 * n_np_y) * n_ed;   % number of equations

n_ee = n_ed * n_en;     % number of element equations

% quadrature rule
n_int_xi  = 3;                     % number of quadrature points in x-direction
n_int_eta = 3;                     % number of quadrature points in y-direction
n_int     = n_int_xi * n_int_eta;  % number of total quadrature points

% =========================================================================
% generate the coordinates of the nodal points
x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);

hh_x = len / (n_el_x * deg_xi);         % mesh size in the x-direction
hh_y = zeros(n_np_x, 1);
for nx = 1 : n_np_x
    hh_y(nx) = (len + heit - hh_x * (nx - 1)) / (n_el_y * deg_eta);        % mesh size in the y-direction
end

for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index for the (nx, ny) node
    x_coor(index) = (nx-1) * hh_x;
    y_coor(index) = (ny-1) * hh_y(nx) + hh_x * (nx - 1) ;
  end
end


% ID and LM arrays are generated based on the BC info
ID = zeros(n_ed, n_np);
counter = 1;
for ny = 1 : n_np_y
    for nx = 1 : n_np_x
        AA = (ny - 1) * n_np_x + nx;
        % Modify ID according to the Dirichlet BC info
        if nx ~= 1 &&  nx ~= n_np_x 
            for ii = 1 : n_ed
                ID(ii, AA) = counter;
                counter = counter + 1;
            end
        end

    end
end


% bilinear quadraliteral element
% setup the IEN array for element with local node numbering as
% a=4 ------- a=3
% |           |
% |           |
% |           |
% |           |
% a=1 ------- a=2
IEN = zeros(n_en, n_el);
for ey = 1 : n_el_y
  for ex = 1 : n_el_x
    ee = (ey-1)*n_el_x + ex;
    IEN(1,ee) = (deg_eta * ey-1)* n_np_x + deg_xi * (ex - 1) + 1;
    IEN(2,ee) = (deg_eta * ey-1)* n_np_x + deg_xi * (ex - 1) + 2;
    IEN(3,ee) =  deg_eta * ey   * n_np_x + deg_xi * (ex - 1) + 2;
    IEN(4,ee) =  deg_eta * ey   * n_np_x + deg_xi * (ex - 1) + 1;
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

% find the node number of Neumann BC
h1_node_number = [];
h2_node_number = [];
for ny = 1 : n_np_y
    for nx = 1 : n_np_x
        AA = (ny - 1) * n_np_x + nx;
        % elements corresponding to the Neumann BC
        if ny == 1 
            h2_node_number(end+1) = AA;
        elseif ny == n_np_y 
            h1_node_number(end+1) = AA;
        end
    end
end

% =========================================================================
% generate the quadrature rule
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);
[xi1D, weight1D] = Gauss(n_int_xi, -1, 1);
% =========================================================================

% Start the assembly procedure
K = spalloc(n_eq, n_eq, 17*n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
   k_ele = zeros(n_ee, n_ee);      % dimension of element stiffness is n_ee x n_ee
   f_ele = zeros(n_ee, 1);         % element force nodes
 
   x_ele = zeros(n_en, 1);         % x-coordinate of nodes in ee th element
   y_ele = zeros(n_en, 1);         % y-coordinate of nodes in ee th element
   for aa = 1 : n_en
     x_ele(aa) = x_coor( IEN(aa,ee) );  % from element node number to global node number 
     y_ele(aa) = y_coor( IEN(aa,ee) );  % A = IEN(a,e)
   end

   % loop over quadrature points   
   for ll = 1 : n_int
     x_l = 0.0; y_l = 0.0;           % coordinate in terms of xi(ll)
     dx_dxi = 0.0; dx_deta = 0.0;
     dy_dxi = 0.0; dy_deta = 0.0;
     % coordinate at current quadrature point
     for aa = 1 : n_en
        x_l = x_l + x_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
        y_l = y_l + y_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
        [Na_xi, Na_eta] = PolyShape_2d_grad(deg_xi, deg_eta, aa, xi(ll), eta(ll));
        dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
        dx_deta = dx_deta + x_ele(aa) * Na_eta;
        dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
        dy_deta = dy_deta + y_ele(aa) * Na_eta;
     end

     detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

     
     const = detJ * weight(ll);
     c1 = lambda_bar + mu;
     c2 = mu;
     c3 = lambda_bar;
     k_iajb = zeros(n_ee, n_ee);

     for aa = 1: n_en
         [Na_xi, Na_eta] = PolyShape_2d_grad(deg_xi, deg_eta, aa, xi(ll), eta(ll));
         Na_x = (Na_xi * dy_deta    - Na_eta * dy_dxi) / detJ;
         Na_y = (Na_xi * (-dx_deta) + Na_eta * dx_dxi)  / detJ;
         Na_grad = [Na_x, Na_y];
         for ii = 1:n_ed
             pp = n_ed*(aa - 1) + ii;
             temp = const * Na_grad(ii);

             f_ele(pp) = f_ele(pp) + weight(ll) * detJ * f{ii}(x_l, y_l) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
             for bb = 1:n_en
                 [Nb_xi, Nb_eta] = PolyShape_2d_grad(deg_xi, deg_eta, bb, xi(ll), eta(ll));
                 Nb_x = (Nb_xi * dy_deta    - Nb_eta * dy_dxi) / detJ;
                 Nb_y = (Nb_xi * (-dx_deta) + Nb_eta * dx_dxi)  / detJ;
                 Nb_grad = [Nb_x, Nb_y];
                 for jj =1:n_ed
                     qq = n_ed *(bb - 1) + jj;
                     k_iajb(pp, qq) = k_iajb(pp,qq) + Nb_grad(jj) * temp;
                 end
             end
         end
     end
     
     for aa = 1:n_en
       for bb = 1 : n_en
           temp = 0;
           for kk = 1 : n_ed
               pp = n_ed * (aa - 1) + kk;
               qq = n_ed * (bb - 1) + kk;
               temp = temp + k_iajb(pp,qq);
           end

           for ii = 1:n_ed
               for jj = 1:n_ed 
                   if ii == jj
                       pp = n_ed * (aa - 1) + ii;
                       qq = n_ed * (bb - 1) + ii;
                       k_ele(pp, qq) = k_ele(pp, qq) + c1 * k_iajb(pp, qq) + c2 * temp;
                   elseif aa == bb
                       pp = n_ed * (aa - 1) + ii;
                       qq = n_ed * (aa - 1) + jj;
                       k_ele(pp, qq) = k_ele(pp, qq) +  c1 * k_iajb(pp, qq);
                   else
                       pp = n_ed * (aa - 1) + ii;
                       qq = n_ed * (bb - 1) + jj;
                       pp_prime = n_ed * (aa - 1) + jj;
                       qq_prime = n_ed * (bb - 1) + ii;
                       k_ele(pp, qq) = k_ele(pp, qq) +  c3 * k_iajb(pp, qq) + c2 * k_iajb(pp_prime, qq_prime);
                   end
               end
           end
       end
     end
   end % end of quadrature loop
   
   

   
   
   % loop over quadrature points for h boundary condtition 
   h_ele = zeros(n_ee, 1);
   for ll = 1:n_int_xi
        x_l_1D_h1 = 0.0;
        x_l_1D_h2 = 0.0;
        dx_dxi_1D_h1 = 0.0;
        dx_dxi_1D_h2 = 0.0;
        for aa = 1: n_en
            if ismember(IEN(aa, ee),h1_node_number)
                
                x_l_1D_h1 = x_l_1D_h1 + x_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi1D(ll), 1);
                [Na_xi, Na_eta] = PolyShape_2d_grad(deg_xi, deg_eta, aa, xi1D(ll), 1);
                dx_dxi_1D_h1  = dx_dxi_1D_h1 + x_ele(aa) * Na_xi;

            elseif ismember(IEN(aa, ee),h2_node_number)
                x_l_1D_h2 = x_l_1D_h2 + x_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi1D(ll), -1);
                [Na_xi, Na_eta] = PolyShape_2d_grad(deg_xi, deg_eta, aa, xi1D(ll), -1);
                dx_dxi_1D_h2  = dx_dxi_1D_h2 + x_ele(aa) * Na_xi;

            end
        end
        for aa = 1:n_en
             for ii = 1 : n_ed
                 pp = n_ed * (aa - 1) + ii;
                 if ismember(IEN(aa, ee),h1_node_number)
                     f_ele(pp) = f_ele(pp) + weight1D(ll) * dx_dxi_1D_h1 * h1{ii}(x_l_1D_h1, len + heit) * PolyShape_2d(deg_xi,deg_eta, aa, xi1D(ll), 1);
                 elseif ismember(IEN(aa, ee),h2_node_number)
                     f_ele(pp) = f_ele(pp) + weight1D(ll) * dx_dxi_1D_h2 * h2{ii}(x_l_1D_h2, x_l_1D_h2) * PolyShape_2d(deg_xi,deg_eta, aa, xi1D(ll), -1);
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
                        F(PP) = F(PP) - k_ele(pp, qq) * g{jj}(x_ele(bb), y_ele(bb));
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
  for aa = 1 : n_en
     x_ele(aa) = x_coor( IEN(aa,ee) );
     y_ele(aa) = y_coor( IEN(aa,ee) );
  end
  for aa = 1:n_en
      for ii = 1: n_ed
          index  = (IEN(aa,ee) - 1) * n_ed + ii;
          pp = n_ed * (aa - 1) + ii;
          PP = LM(pp, ee);           % equation number 
          if PP > 0
              disp(index) = d_temp(PP);
          else
              disp(index) = g{ii}(x_ele(aa), y_ele(aa));
          end
      end
   end
end

% Extract the displacement component
ux_h = zeros(n_np * n_ed / 2, 1);
% Extract the angle component
uy_h = zeros(n_np * n_ed / 2, 1);
for ii = 1 : n_np * n_ed
    if mod(ii, 2) == 1
        ux_h(floor(ii / 2) + 1) = disp(ii);
    else
        uy_h(floor(ii / 2)) = disp(ii);
    end
end



% plot mesh
figure;
for ee = 1 : n_el
    XX = [x_coor(IEN(1, ee)), x_coor(IEN(2, ee)), x_coor(IEN(3, ee)), x_coor(IEN(4, ee)), x_coor(IEN(1, ee))];
    YY = [y_coor(IEN(1, ee)), y_coor(IEN(2, ee)), y_coor(IEN(3, ee)), y_coor(IEN(4, ee)), y_coor(IEN(1, ee))];
    plot (XX, YY,'LineWidth', 1);
    hold on;
end
% plot displacement
xnew_coor = zeros(n_np, 1);
ynew_coor = zeros(n_np, 1);
for ii = 1 :  n_np
    xnew_coor(ii) = x_coor(ii) + ux_h(ii);
    ynew_coor(ii) = y_coor(ii) + uy_h(ii);
end
for ee = 1 : n_el
    XXnew = [xnew_coor(IEN(1, ee)), xnew_coor(IEN(2, ee)), xnew_coor(IEN(3, ee)), xnew_coor(IEN(4, ee)), xnew_coor(IEN(1, ee))];
    YYnew = [ynew_coor(IEN(1, ee)), ynew_coor(IEN(2, ee)), ynew_coor(IEN(3, ee)), ynew_coor(IEN(4, ee)), ynew_coor(IEN(1, ee))];
    plot (XXnew, YYnew, "r-", 'LineWidth', 1);
    hold on;
end
title("Undeformed v.s. Deformed");





% plot the solution ux
figure;
% FEM solution
X_h = reshape(x_coor, n_np_x, n_np_y);
Y_h = reshape(y_coor, n_np_x, n_np_y);
Z_h = reshape(ux_h, n_np_x, n_np_y);
surf1 = surf(X_h, Y_h, Z_h);
set(surf1, 'FaceColor', 'red');
hold on;
% exact solution
xexact_coor = zeros(101 * 101, 1);
yexact_coor = zeros(101 * 101, 1);
hh_x = len / 100;         % mesh size in the x-direction
hh_y = zeros(100, 1);
for nx = 1 : 101
    hh_y(nx) = (len + heit - hh_x * (nx - 1)) / 100;        % mesh size in the y-direction
end

for ny = 1 : 101
  for nx = 1 : 101
    index = (ny-1)*101 + nx; % nodal index for the (nx, ny) node
    xexact_coor(index) = (nx-1) * hh_x;
    yexact_coor(index) = (ny-1) * hh_y(nx) + hh_x * (nx - 1) ;
  end
end
X = reshape(xexact_coor, 101, 101);
Y = reshape(yexact_coor, 101, 101);
for ii = 1:101
    for jj = 1:101
        Z(ii,jj) = exact_ux(X(ii,jj), Y(ii,jj));
    end
end
surf2 = surf(X, Y, Z);
set(surf2, 'FaceColor', 'blue');
xlabel("X");
ylabel("Y");
zlabel("Displacement");
legend('有限元解-x-位移', '精确解-x-位移','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% plot the solution uy
figure;
% FEM solution
X_h = reshape(x_coor, n_np_x, n_np_y);
Y_h = reshape(y_coor, n_np_x, n_np_y);
Z_h = reshape(uy_h, n_np_x, n_np_y);
surf1 = surf(X_h, Y_h, Z_h);
set(surf1, 'FaceColor', 'red');
hold on;

for ii = 1:101
    for jj = 1:101
        Z(ii,jj) = exact_uy(X(ii,jj), Y(ii,jj));
    end
end
surf2 = surf(X, Y, Z);
set(surf2, 'FaceColor', 'blue');
xlabel("X");
ylabel("Y");
zlabel("Displacement");
legend('有限元解-y-位移', '精确解-y-位移','Location', 'Best', 'FontSize', 14, 'Box', 'on');




% postprocess the solution by calculating the error measured in L2 norm
nqp = 10; % we need more points 
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

errorL2 = 0.0; bottomL2 = 0.0;
errorH1 = 0.0; bottomH1 = 0.0;
for ee = 1 : n_el
  x_ele = x_coor( IEN(1:n_en, ee) );
  y_ele = y_coor( IEN(1:n_en, ee) );
  ux_ele = ux_h(   IEN(1:n_en, ee) );
  uy_ele = uy_h(   IEN(1:n_en, ee) );

  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0; ux_l = 0.0; uy_l = 0.0;
    ux_l_xi = 0.0; ux_l_eta = 0.0;uy_l_xi = 0.0; uy_l_eta = 0.0;
    dx_dxi = 0.0; dy_dxi = 0.0; dx_deta = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
      ux_l = ux_l + ux_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
      uy_l = uy_l + uy_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = PolyShape_2d_grad(deg_xi, deg_eta, aa, xi(ll), eta(ll));
      ux_l_xi  = ux_l_xi  + ux_ele(aa) * Na_xi;
      uy_l_xi  = uy_l_xi  + uy_ele(aa) * Na_xi;
      ux_l_eta = ux_l_eta + ux_ele(aa) * Na_eta;
      uy_l_eta = uy_l_eta + uy_ele(aa) * Na_eta;
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

    ux_l_x = (ux_l_xi * dy_deta - ux_l_eta * dy_dxi) / detJ;
    ux_l_y = (ux_l_xi * (-dx_deta) + ux_l_eta * dx_dxi) / detJ;
    uy_l_x = (uy_l_xi * dy_deta - uy_l_eta * dy_dxi) / detJ;
    uy_l_y = (uy_l_xi * (-dx_deta) + uy_l_eta * dx_dxi) / detJ;

    errorL2 = errorL2 + weight(ll) * detJ * ( (ux_l - exact_ux(x_l, y_l))^2 + ...
        (uy_l - exact_uy(x_l, y_l))^2 );
    errorH1 = errorH1 + weight(ll) * detJ *...
      (( ux_l_x- exact_ux_x(x_l,y_l))^2 + ( ux_l_y - exact_ux_y(x_l,y_l))^2 + ... 
       ( uy_l_x- exact_uy_x(x_l,y_l))^2 + ( uy_l_y - exact_uy_y(x_l,y_l))^2 );
    bottomL2 = bottomL2 + weight(ll) * detJ * (exact_ux(x_l, y_l)^2 + exact_uy(x_l, y_l)^2);
    bottomH1 = bottomH1 + weight(ll) * detJ * (exact_ux_x(x_l,y_l)^2 + exact_ux_y(x_l,y_l)^2 ...
        + exact_uy_x(x_l,y_l)^2 + exact_uy_y(x_l,y_l)^2);
  end
end

errorL2 = sqrt(errorL2) / sqrt(bottomL2);
errorH1 = sqrt(errorH1) / sqrt(bottomH1);

% EOF