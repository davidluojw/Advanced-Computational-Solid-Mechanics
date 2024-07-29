clear all; clc;

kappa = 1.0; % isotropic homogeneous heat conductivity
rho   = 1.0; % density
cap   = 1.0; % heat capacity

% manufactured solution and source term
exact   = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

%f = @(x,y) -2*x*(x-1)-2*y*(y-1);

alpha = 0.5;

T_final = 10.0;
dt      = 0.1;
t       = 0 : dt : T_final; % time sub-interval
NN      = T_final / dt;

% quadrature rule
n_int_xi  = 3;              % number of quadrature points in xi-direction
n_int_eta = 3;              % number of quadrature points in eta-direction
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% FEM mesh settings
n_en = 4;                   % 4-node quadrilateral element

n_el_x = 200;               % number of element in x-direction
n_el_y = 200;               % number of element in y-direction
n_el   = n_el_x * n_el_y;   % total number of element in 2D domain

n_np_x = n_el_x + 1;        % number of node points in x-direction
n_np_y = n_el_y + 1;        % number of node points in y-direction
n_np   = n_np_x * n_np_y;   % total number of node points in 2D domain

% generate the coordinates of the nodal points
x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);

hh_x = 1 / n_el_x;          % mesh size in the x-direction
hh_y = 1 / n_el_y;          % mesh size in the y-direction

for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index for the (nx, ny) node
    x_coor(index) = (nx-1) * hh_x;
    y_coor(index) = (ny-1) * hh_y;
  end
end

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
    IEN(1,ee) = (ey-1)* n_np_x + ex;
    IEN(2,ee) = (ey-1)* n_np_x + ex + 1;
    IEN(3,ee) =  ey   * n_np_x + ex + 1;
    IEN(4,ee) =  ey   * n_np_x + ex;
  end
end

% ID array
ID = zeros(n_np, 1);
counter = 1;
for ny = 2 : n_np_y-1
  for nx = 2 : n_np_x-1
    ID( (ny-1)*n_np_x + nx ) = counter;
    counter = counter + 1;
  end
end

n_eq = n_np - n_np_x * 2 - n_np_y * 2 + 4;

LM = ID(IEN);

% Start the assembly procedure
K = spalloc(n_eq, n_eq, 9*n_eq);
M = K;
F = zeros(n_eq, 1);

for ee = 1 : n_el
   k_ele = zeros(n_en, n_en);
   m_ele = zeros(n_en, n_en);
   f_ele = zeros(n_en, 1);

   x_ele = x_coor( IEN(1:n_en, ee) );
   y_ele = y_coor( IEN(1:n_en, ee) );

   % loop over quadrature points   
   for ll = 1 : n_int
     % we need the geometric mapping at each quadrature points
     x_l = 0.0; y_l = 0.0;
     dx_dxi = 0.0; dx_deta = 0.0;
     dy_dxi = 0.0; dy_deta = 0.0;
     for aa = 1 : n_en
        x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
        y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
        [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
        dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
        dx_deta = dx_deta + x_ele(aa) * Na_eta;
        dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
        dy_deta = dy_deta + y_ele(aa) * Na_eta;
     end

     detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

     % we loop over a and b to assemble the element stiffness matrix and load
     % vector
     for aa = 1 : n_en
       Na = Quad(aa, xi(ll), eta(ll));
       [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));

       % See page 147 of Sec. 3.9 of the textbook for the shape function
       % routine details
       Na_x = (Na_xi * dy_deta    - Na_eta * dy_dxi) / detJ;
       Na_y = (Na_xi * (-dx_deta) + Na_eta * dx_dxi)  / detJ;

       f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Quad(aa, xi(ll), eta(ll));
       for bb = 1 : n_en
         Nb = Quad(bb, xi(ll), eta(ll));
         [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
         Nb_x = (Nb_xi * dy_deta    - Nb_eta * dy_dxi) / detJ;
         Nb_y = (Nb_xi * (-dx_deta) + Nb_eta * dx_dxi)  / detJ;

         m_ele(aa,bb) = m_ele(aa,bb) + weight(ll) * detJ * rho * cap * Na * Nb;
         k_ele(aa,bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * ( Na_x * Nb_x + Na_y * Nb_y);

       end % end of bb-loop
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
           M(PP, QQ) = M(PP, QQ) + m_ele(aa, bb);
           K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
         else
           % do something for non-zero g boundary condition
         end
       end
     end
   end
end % end of element loop

% set initial conditions
dn = zeros(n_eq, 1);   % initial temperature is zero
vn = M \ (F - K * dn); % solve M matrix to determine vn at initial time
disp = zeros(n_np, 1);

% insert the solution vector back with the g-data
for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    disp(ii) = dn(index);
  end
end
save("HEAT"+int2str(1000000), "disp", "n_el_x", "n_el_y");

LEFT = M + alpha * dt * K;

for n = 1 : NN
  % prediction
  tilde_d = dn + (1-alpha) * dt * vn;
  
  % correction
  RIGHT = F - K * tilde_d;
  vn = LEFT \ RIGHT;
  dn = tilde_d + alpha * dt * vn;
 
  disp = zeros(n_np, 1);

  % insert the solution vector back with the g-data
  for ii = 1 : n_np
    index = ID(ii);
    if index > 0
      disp(ii) = dn(index);
    end
  end
  
  % save the solution to file
  save("HEAT"+int2str(1000000+n), "disp", "n_el_x", "n_el_y");
end

% EOF