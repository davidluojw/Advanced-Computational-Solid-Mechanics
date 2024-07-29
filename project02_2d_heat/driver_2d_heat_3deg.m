clear all; clc;

% =========================================================================
% Problem definition
% exact solution
% manufactured solution and source term
G = @(x, y) sin((x + y) * 2 * pi);
G_x = @(x, y) 2 * pi * cos((x + y) * 2 * pi);
G_y = @(x, y) 2 * pi * cos((x + y) * 2 * pi);
G_xx = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);
G_yy = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);
G_xy = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);

exact   = @(x,y) x*(1-x)*y*(1-y) + 0.1*G(x, y);
exact_x = @(x,y) (1-2*x)*y*(1-y) + 0.1*G_x(x, y); 
exact_y = @(x,y) x*(1-x)*(1-2*y) + 0.1*G_y(x, y);
exact_xx = @(x,y) -2*y*(1-y) + 0.1*G_xx(x, y); 
exact_yy = @(x,y) -2*x*(1-x) + 0.1*G_yy(x, y);
exact_xy = @(x,y) (1-2*x)*(1-2*y) + 0.1*G_xy(x, y); 

kappa = eye(2); % isotropic homogeneous heat conductivity

% force
f = @(x,y) -kappa(1,1) * exact_xx(x,y) -2 * kappa(1,2) * exact_xy(x,y) - kappa(2,2) * exact_yy(x,y);



% Dirichlet BC
g = @(x, y) 0.1*G(x, y);


% Neumann BC
% h_0 = @(x, y) -kappa(2,1) * 0.1 * G(x,y) - kappa(2,2) * exact_y(x, y);
% h_1 = @(x, y) kappa(2,1) * 0.1 * G(x,y) + kappa(2,2) * exact_y(x, y);
h_0 = @(x,y) 0;
h_1 = @(x,y) 0;


% =========================================================================

% =========================================================================
% Generate the mesh
% FEM mesh settings
n_sd = 2;                 % space dimension

n_el_x = 7;               % number of element in x-direction
n_el_y = 6;               % number of element in y-direction
n_el   = n_el_x * n_el_y;   % total number of element in 2D domain

n_en = 16; % number of  element nodes
n_ed = 1; % number of element degrees of freedom (per node)

deg_xi = 3;   % degree of lagrange polynomial in xi-direction
deg_eta = 3;  % degree of lagrange polynomial in eta-direction

n_np_x = n_el_x * deg_xi  + 1;      % number of node points in x-direction
n_np_y = n_el_y * deg_eta + 1;      % number of node points in y-direction
n_np   = n_np_x * n_np_y; % total number of node points in 2D domain

n_eq = n_np - n_np_y * 2;   % number of equations

n_ee = n_ed * n_en;     % number of element equations

% quadrature rule
n_int_xi  = 3;                     % number of quadrature points in x-direction
n_int_eta = 3;                     % number of quadrature points in y-direction
n_int     = n_int_xi * n_int_eta;  % number of total quadrature points

% =========================================================================
% generate the coordinates of the nodal points
x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);

hh_x = 1 / (n_el_x * deg_xi);         % mesh size in the x-direction
hh_y = 1 / (n_el_y * deg_eta);        % mesh size in the y-direction

for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index for the (nx, ny) node
    x_coor(index) = (nx-1) * hh_x;
    y_coor(index) = (ny-1) * hh_y;
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



% biquadratic quadraliteral element
% setup the IEN array for element with local node numbering as
% a=4 ---11-----7-----a=3
% |      |      |     |
% |      |      |     |
% 8------16-----15----10
% |      |      |     |
% |      |      |     |
% 12-----13-----14----6
% |      |      |     |
% |      |      |     |
% a=1 ---5------9-----a=2
IEN = zeros(n_en, n_el);
for ey = 1 :  n_el_y
  for ex = 1 : n_el_x
    ee = (ey-1)*n_el_x + ex;
    IEN(1,ee) = (deg_eta * ey-3)* n_np_x + deg_xi * (ex - 1) + 1;
    IEN(2,ee) = (deg_eta * ey-3)* n_np_x + deg_xi * (ex - 1) + 4;
    IEN(3,ee) =  deg_eta * ey   * n_np_x + deg_xi * (ex - 1) + 4;
    IEN(4,ee) =  deg_eta * ey   * n_np_x + deg_xi * (ex - 1) + 1;
    IEN(5,ee) = (deg_eta * ey-3)* n_np_x + deg_xi * (ex - 1) + 2;
    IEN(6,ee) = (deg_eta * ey-2)* n_np_x + deg_xi * (ex - 1) + 4;
    IEN(7,ee) =  deg_eta * ey   * n_np_x + deg_xi * (ex - 1) + 3;
    IEN(8,ee) = (deg_eta * ey-1)* n_np_x + deg_xi * (ex - 1) + 1;
    IEN(9,ee) = (deg_eta * ey-3)* n_np_x + deg_xi * (ex - 1) + 3;
    IEN(10,ee)= (deg_eta * ey-1)* n_np_x + deg_xi * (ex - 1) + 4;
    IEN(11,ee)=  deg_eta * ey   * n_np_x + deg_xi * (ex - 1) + 2;
    IEN(12,ee)= (deg_eta * ey-2)* n_np_x + deg_xi * (ex - 1) + 1;
    IEN(13,ee)= (deg_eta * ey-2)* n_np_x + deg_xi * (ex - 1) + 2;
    IEN(14,ee)= (deg_eta * ey-2)* n_np_x + deg_xi * (ex - 1) + 3;
    IEN(15,ee)= (deg_eta * ey-1)* n_np_x + deg_xi * (ex - 1) + 3;
    IEN(16,ee)= (deg_eta * ey-1)* n_np_x + deg_xi * (ex - 1) + 2;
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

h0_node_number = [];
h1_node_number = [];
for ny = 1 : n_np_y
    for nx = 1 : n_np_x
        AA = (ny - 1) * n_np_x + nx;

        % elements corresponding to the Neumann BC
        if ny == 1
            h0_node_number(end+1) = AA;
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
K = spalloc(n_eq, n_eq, 9*n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
   k_ele = zeros(n_ee, n_ee);      % dimension of element stiffness is n_ee x n_ee
   f_ele = zeros(n_ee, 1);         % element force nodes
 
   x_ele = zeros(n_en, 1);         % coordinate of nodes in ee th element
   y_ele = x_ele;
   for aa = 1 : n_en
     x_ele(aa) = x_coor( IEN(aa,ee) );  % from element node number to global node number 
     y_ele(aa) = y_coor( IEN(aa,ee) );  % A = IEN(a,e)
   end

   % loop over quadrature points   
   for ll = 1 : n_int
     x_l = 0.0; y_l = 0.0;           % coordinate in terms of xi(ll)
     dx_dxi = 0.0; dx_deta = 0.0;
     dy_dxi = 0.0; dy_deta = 0.0;
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
     
     for aa = 1 : n_en
       [Na_xi, Na_eta] = PolyShape_2d_grad(deg_xi, deg_eta, aa, xi(ll), eta(ll));
       Na_x = (Na_xi * dy_deta    - Na_eta * dy_dxi) / detJ;
       Na_y = (Na_xi * (-dx_deta) + Na_eta * dx_dxi)  / detJ;

       f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
       for bb = 1 : n_en
         [Nb_xi, Nb_eta] = PolyShape_2d_grad(deg_xi, deg_eta, bb, xi(ll), eta(ll));
         Nb_x = (Nb_xi * dy_deta    - Nb_eta * dy_dxi) / detJ;
         Nb_y = (Nb_xi * (-dx_deta) + Nb_eta * dx_dxi)  / detJ;

         k_ele(aa,bb) = k_ele(aa,bb) + weight(ll) * detJ * [Na_x, Na_y] * kappa * [Nb_x; Nb_y];
       end % end of bb-loop
     end % end of aa-loop
   end % end of quadrature loop
   

   % loop over quadrature points for h boundary condtition 
   for ll = 1:n_int_xi
        x_l_1D_n1 = 0.0; 
        x_l_1D_p1 = 0.0;
        dx_dxi_1D_n1 = 0.0;
        dx_dxi_1D_p1 = 0.0;
        for aa = 1: n_en
            if ismember(IEN(aa, ee),h0_node_number)
                x_l_1D_n1 = x_l_1D_n1 + x_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi1D(ll), -1);
                [Na_xi, Na_eta] = PolyShape_2d_grad(deg_xi, deg_eta, aa, xi1D(ll), -1);
                dx_dxi_1D_n1  = dx_dxi_1D_n1 + x_ele(aa) * Na_xi;
            elseif ismember(IEN(aa, ee),h1_node_number)
                x_l_1D_p1 = x_l_1D_p1 + x_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi1D(ll), 1);
                [Na_xi, Na_eta] = PolyShape_2d_grad(deg_xi, deg_eta, aa, xi1D(ll), 1);
                dx_dxi_1D_p1  = dx_dxi_1D_p1 + x_ele(aa) * Na_xi;
            end
        end
        for aa = 1:n_en
           if ismember(IEN(aa, ee),h0_node_number)
               f_ele(aa) = f_ele(aa) + weight1D(ll) * dx_dxi_1D_n1 * h_0(x_l_1D_n1, 0) * PolyShape_2d(deg_xi,deg_eta, aa, xi1D(ll), -1);
           elseif ismember(IEN(aa, ee),h1_node_number)
               f_ele(aa) = f_ele(aa) + weight1D(ll) * dx_dxi_1D_p1 * h_1(x_l_1D_p1, 1) * PolyShape_2d(deg_xi,deg_eta, aa, xi1D(ll), 1);
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
           % g: {x = 0, y in (0,1)}+{x = 1, y in (0,1)}
           F(PP) = F(PP) - k_ele(aa, bb) * g(x_ele(bb), y_ele(bb));
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
  for aa = 1 : n_en
     x_ele(aa) = x_coor( IEN(aa,ee) );
     y_ele(aa) = y_coor( IEN(aa,ee) );
  end
  
  for aa = 1:n_en
    index = LM(aa, ee);           % equation number 
    AA = IEN(aa,ee);              % global node number 
    if index > 0
        disp(AA) = d_temp(index);
    else
        disp(AA) = g(x_ele(aa), y_ele(aa));
    end
  end
end

% plot the solution
figure;
[X_h, Y_h] = meshgrid( 0:hh_x:1, 0:hh_y:1 );
Z_h = reshape(disp, n_np_x, n_np_y);
surf(X_h, Y_h, Z_h');
hold on;
[X, Y] = meshgrid( 0:0.01:1, 0:0.01:1 );
for ii = 1:101
    for jj = 1:101
        Z(ii,jj) = exact(X(ii,jj), Y(ii,jj));
    end
end
surf(X, Y, Z);
xlabel("X");
ylabel("Y");
zlabel("Temperature");
legend('有限元解-温度', '精确解-温度','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% plot the solution
figure;
[X_h, Y_h] = meshgrid( 0:hh_x:1, 0:hh_y:1 );
Z_h = reshape(disp, n_np_x, n_np_y);
surf(X_h, Y_h, Z_h');
xlabel("X");
ylabel("Y");
zlabel("Temperature");
legend('有限元解-温度','Location', 'Best', 'FontSize', 14, 'Box', 'on');
figure;
[X, Y] = meshgrid( 0:0.01:1, 0:0.01:1 );
for ii = 1:101
    for jj = 1:101
        Z(ii,jj) = exact(X(ii,jj), Y(ii,jj));
    end
end
surf(X, Y, Z);
xlabel("X");
ylabel("Y");
zlabel("Temperature");
legend('精确解-温度','Location', 'Best', 'FontSize', 14, 'Box', 'on');


% postprocess the solution by calculating the error measured in L2 norm
nqp = 10; % we need more points 
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

errorL2 = 0.0; bottomL2 = 0.0;
errorH1 = 0.0; bottomH1 = 0.0;
for ee = 1 : n_el
  x_ele = x_coor( IEN(1:n_en, ee) );
  y_ele = y_coor( IEN(1:n_en, ee) );
  u_ele = disp(   IEN(1:n_en, ee) );

  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0; u_l = 0.0;
    u_l_xi = 0.0; u_l_eta = 0.0;
    dx_dxi = 0.0; dy_dxi = 0.0; dx_deta = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
      u_l = u_l + u_ele(aa) * PolyShape_2d(deg_xi, deg_eta, aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = PolyShape_2d_grad(deg_xi, deg_eta, aa, xi(ll), eta(ll));
      u_l_xi  = u_l_xi  + u_ele(aa) * Na_xi;
      u_l_eta = u_l_eta + u_ele(aa) * Na_eta;
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

    u_l_x = (u_l_xi * dy_deta - u_l_eta * dy_dxi) / detJ;
    u_l_y = (u_l_xi * (-dx_deta) + u_l_eta * dx_dxi) / detJ;

    errorL2 = errorL2 + weight(ll) * detJ * (u_l - exact(x_l, y_l))^2;
    errorH1 = errorH1 + weight(ll) * detJ *...
      (( u_l_x- exact_x(x_l,y_l))^2 + ( u_l_y - exact_y(x_l,y_l))^2);
    bottomL2 = bottomL2 + weight(ll) * detJ * exact(x_l, y_l)^2;
    bottomH1 = bottomH1 + weight(ll) * detJ * (exact_x(x_l,y_l)^2 + exact_y(x_l,y_l)^2);
  end
end

errorL2 = sqrt(errorL2) / sqrt(bottomL2);
errorH1 = sqrt(errorH1) / sqrt(bottomH1);

% EOF