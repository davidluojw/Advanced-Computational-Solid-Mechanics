% This is the nonlinear solver for heat equation with material nonlinearity
% clean the memory and the screen
clear all; clc;

% domain
omega_l = 0.0;
omega_r = 1.0;

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
h = @(x) -1.0 * fun_kappa( exact(0) );
g = @(x) exact(1);
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
uh = [ zeros(n_eq,1); g(omega_r) ];

counter = 0;
nmax    = 200;
error   = 1.0;

% Allocate an empty stiffness matrix and load vector
K = sparse(n_eq, n_eq);
F = zeros(n_eq, 1);
  
% Assembly the siffness matrix and load vector
for ee = 1 : nElem
    % Allocate zero element stiffness matrix and element load vector
    k_ele = zeros(n_en, n_en);
    f_ele = zeros(n_en, 1);
    
    x_ele = zeros(n_en, 1);
    d_ele = zeros(n_en, 1);
    for aa = 1 : n_en
        x_ele(aa) = x_coor( IEN(aa, ee) );
        u_ele(aa) = uh( IEN(aa, ee) );
    end
    
    for qua = 1 : nqp
        % geometrical mapping
        x_qua    = 0.0;
        dx_dxi   = 0.0;
        u_qua    = 0.0;
        u_xi     = 0.0;
        for aa = 1 : n_en
            x_qua    = x_qua  + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            dx_dxi   = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
            u_qua    = u_qua  + u_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            u_xi     = u_xi   + u_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
        end
        dxi_dx = 1.0 / dx_dxi;
          
        kappa = fun_kappa( u_qua );
        dkappa = fun_dkappa( u_qua );
          
        for aa = 1 : n_en
            Na    = PolyBasis(pp, aa, 0, qp(qua));
            Na_xi = PolyBasis(pp, aa, 1, qp(qua));
            f_ele(aa) = f_ele(aa) + wq(qua) * Na * f(x_qua) * dx_dxi;
            f_ele(aa) = f_ele(aa) - wq(qua) * Na_xi * kappa * u_xi * dxi_dx;
            
            for bb = 1 : n_en
            Nb    = PolyBasis(pp, bb, 0, qp(qua));
            Nb_xi = PolyBasis(pp, bb, 1, qp(qua));
            k_ele(aa,bb) = k_ele(aa,bb) + wq(qua) * Na_xi * kappa * Nb_xi * dxi_dx;
            k_ele(aa,bb) = k_ele(aa,bb) + wq(qua) * Na_xi * dkappa * Nb * u_xi * dxi_dx;
            end
        end
    end
    % end of the quadrature loop
    
    % distribute the entries to the global stiffness matrix and global load vector
    for aa = 1 : n_en
        LM_a = ID( IEN(aa, ee) );
        if LM_a > 0
            F(LM_a) = F(LM_a) + f_ele(aa);
            for bb = 1 : n_en
                LM_b = ID( IEN(bb, ee) );
                if LM_b > 0
                  K(LM_a, LM_b) = K(LM_a, LM_b) + k_ele(aa, bb);
                else
                  % x_qua = x_coor( IEN(bb,ee) ); % obtain the Dirichlet node's physical coordinates
                  % g_qua = g( x_qua ); % Obtain the boundary data at this point
                  % F( LM_a ) = F( LM_a ) - k_ele(aa, bb) * g( x_ele(bb) );
                  % 这里不应该有F的组装，因为这个方程解的是一个增量，没有Dirichlet的影响
                end
            end
        end
    end
    
    % Modify the load vector by the Natural BC
    % Note: for multi-dimensional cases, one needs to perform line or
    % surface integration for the natural BC.
    if ee == 1
        F( ID(IEN(1, ee)) ) = F( ID(IEN(1, ee)) ) + h( x_coor(IEN(1,ee)));
    end
end

while counter < nmax && error > 1.0e-8

    % Solve the stiffness matrix problem
    incremental = K \ F;
    uh = [ uh(1:end-1) + incremental; g(omega_r) ];

    G0 = incremental' * F;
    
    % Assembly the siffness matrix and load vector
    for ee = 1 : nElem
        % Allocate zero element stiffness matrix and element load vector
        k_ele = zeros(n_en, n_en);
        f_ele = zeros(n_en, 1);
        
        x_ele = zeros(n_en, 1);
        d_ele = zeros(n_en, 1);
        for aa = 1 : n_en
            x_ele(aa) = x_coor( IEN(aa, ee) );
            u_ele(aa) = uh( IEN(aa, ee) );
        end
        
        for qua = 1 : nqp
            % geometrical mapping
            x_qua    = 0.0;
            dx_dxi   = 0.0;
            u_qua    = 0.0;
            u_xi     = 0.0;
            for aa = 1 : n_en
                x_qua    = x_qua  + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
                dx_dxi   = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
                u_qua    = u_qua  + u_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
                u_xi     = u_xi   + u_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
            end
            dxi_dx = 1.0 / dx_dxi;
              
            kappa = fun_kappa( u_qua );
            dkappa = fun_dkappa( u_qua );
              
            for aa = 1 : n_en
                Na    = PolyBasis(pp, aa, 0, qp(qua));
                Na_xi = PolyBasis(pp, aa, 1, qp(qua));
                f_ele(aa) = f_ele(aa) + wq(qua) * Na * f(x_qua) * dx_dxi;
                f_ele(aa) = f_ele(aa) - wq(qua) * Na_xi * kappa * u_xi * dxi_dx;
             
            end
        end
        % end of the quadrature loop
        
        % distribute the entries to the global stiffness matrix and global load vector
        for aa = 1 : n_en
            LM_a = ID( IEN(aa, ee) );
            if LM_a > 0
                F(LM_a) = F(LM_a) + f_ele(aa);
            end
        end
        
        % Modify the load vector by the Natural BC
        % Note: for multi-dimensional cases, one needs to perform line or
        % surface integration for the natural BC.
        if ee == 1
            F( ID(IEN(1, ee)) ) = F( ID(IEN(1, ee)) ) + h( x_coor(IEN(1,ee)));
        end
    end

    G = increment' * F;

    s = 1;

    % Line Search
    % find a s such that abs(G) <= stol * abs(G0)
    if abs(G) > stol * abs(G0)
        linmax = 10;
        smax = 16;
        sb = 0;
        sa = s;
        Gb = G0;
        Ga = G;

        % temporary variables
        d_temp = d;
        R_temp = zeros(2,1);

        % find bracket on zero
        % if Ga * Gb > 0, [sb,sa] no zero point
        % then enlarge the bracket to [sa, 2sa]
        % compute the corresponding R_temp and Gb, Ga 
        % in new bracket [sa,2sa]
        while (Ga * Gb > 0 && sa < smax)
            sb = sa;
            sa = 2*sa;
            Gb = Ga;
            d_temp = d_temp + sa * Deltad;
            R_temp(1) = R1(d_temp(1),d_temp(2),n);
            R_temp(2) = R2(d_temp(1),d_temp(2),n);
            Ga = Deltad' * R_temp;
        end

        step = sa; 
        G = Ga;
        lin_step = 0;
        % now we have already found the bracket [sb,sa]
        % which contains the zero point. Ga * Gb < 0

        % Illinois algorithm to find zero
        % while still abs(G) > stol * abs(G0), criteria is not satisfied
        while (lin_step <= linmax && Ga * Gb < 0 && ...
                abs(G) > stol * abs(G0) ) %|| abs(sb-sa) > stol * 0.5 * (sb+sa)))
            step = sa - Ga * (sa - sb) / (Ga - Gb);   % linear interpolation
            d_temp = d_temp + step * Deltad;
            R_temp(1) = R1(d_temp(1),d_temp(2),n);
            R_temp(2) = R2(d_temp(1),d_temp(2),n);
            G = Deltad' * R_temp;

            if G * Ga > 0
                Gb = 0.5 * Gb;
            else
                sb = sa;
                Gb = Ga;
            end
            sa = step;
            Ga = G;
            lin_step = lin_step + 1;
        end
        s = step;
        % update the new updates
        d = d - Deltad;
        d = d + s * Deltad;

        R(1) = R1(d(1),d(2),n);
        R(2) = R2(d(1),d(2),n);
        R_norm = sqrt(R(1) * R(1) + R(2) * R(2));
    end

  
    error = norm(F);
    counter = counter + 1;
end

% plot the solution
X_h = 0: hh/pp :1;
Y_h = uh;
% subplot(2,1,1)
plot(X_h, Y_h,'b-', 'LineWidth', 2);
hold on;
% subplot(2,1,2)
X = 0:0.01:1;
Y = exact(X);
plot(X, Y,'r-', 'LineWidth', 2);
xlabel("X");
ylabel("Temperature");
legend('FEM', 'EXACT','Location', 'Best', 'FontSize', 14, 'Box', 'on');

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

error = top / bot;

top_u = sqrt(top_u);
bot_u = sqrt(bot_u);
error_u = top_u / bot_u;

% EOF