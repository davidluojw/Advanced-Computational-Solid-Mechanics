clear all; clc;

load("FEM_solution.mat");

n_el = n_el_x * n_el_y;
n_en = 4;

n_int_x = 10;
n_int_y = 10;
n_int   = n_int_x * n_int_y;

n_np_x = n_el_x + 1;
n_np_y = n_el_y + 1;
n_np   = n_np_x * n_np_y;

n_eq = n_np - n_np_x - n_np_x - n_np_y - n_np_y + 4;

hh_x = 1.0 / n_el_x;
hh_y = 1.0 / n_el_y;

x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);

for nx = 1 : n_np_x
    for ny = 1 : n_np_y
        x_coor( (ny-1)*n_np_x + nx ) = (nx-1)*hh_x;
        y_coor( (ny-1)*n_np_x + nx ) = (ny-1)*hh_y;
    end
end

IEN = zeros(n_en, n_el);
for ex = 1 : n_el_x
    for ey = 1 : n_el_y
        ee = (ey-1) * n_el_x + ex;
        IEN(1, ee) = (ey-1)*n_np_x + ex;
        IEN(2, ee) = (ey-1)*n_np_x + ex + 1;
        IEN(3, ee) = ey * n_np_x + ex + 1;
        IEN(4, ee) = ey * n_np_x + ex;
    end
end

[ xi, eta, weight ] = Gauss2D(n_int_x, n_int_y);

errorL2 = 0.0;
bottomL2 = 0.0;

errorH1 = 0.0;
bottomH1 = 0.0;

for ee = 1 : n_el
    x_ele = x_coor( IEN(1:n_en, ee) );
    y_ele = y_coor( IEN(1:n_en, ee) );
    u_ele = disp( IEN(1:n_en, ee) );

    for ll = 1 : n_int
        x_l = 0.0;
        y_l = 0.0;
        u_l = 0.0;
        u_l_xi = 0.0;
        u_l_eta = 0.0;
        dx_dxi  = 0.0;
        dx_deta = 0.0;
        dy_dxi  = 0.0;
        dy_deta = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            u_l = u_l + u_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
            u_l_xi  = u_l_xi + u_ele(aa) * Na_xi;
            u_l_eta = u_l_eta + u_ele(aa) * Na_eta;
        end
        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

        u_l_x = (u_l_xi * dy_deta - u_l_eta * dy_dxi) / detJ;
        u_l_y = (-u_l_xi * dx_deta + u_l_eta * dx_dxi) / detJ;

        errorL2  = errorL2  + weight(ll) * detJ * ( exact(x_l, y_l) - u_l )^2;
        bottomL2 = bottomL2 + weight(ll) * detJ * exact(x_l, y_l)^2;

        errorH1  = errorH1  + weight(ll) * detJ * ( ( exact_x(x_l, y_l) - u_l_x)^2 + ( exact_y(x_l, y_l) - u_l_y)^2 );
        bottomH1 = bottomH1 + weight(ll) * detJ * ( exact_x(x_l, y_l)^2 + exact_y(x_l, y_l)^2 );
    end
end

L2error = sqrt( errorL2 ) / sqrt(bottomL2);
H1error = sqrt( errorH1 ) / sqrt(bottomH1);

% EOF