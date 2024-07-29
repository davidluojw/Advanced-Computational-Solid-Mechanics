function [val_xi, val_eta, val_zeta]= PolyShape_3d_grad(deg_xi, deg_eta, deg_zeta, a, xi, eta, zeta)

if deg_xi == deg_eta && deg_eta == deg_zeta
    degree = deg_xi;
else
    error('Error: the three degree of xi, eta and zeta must be identical.');
end

switch degree
    % trilinear hexahedral 8-nodes
    case 1
        if a == 1
            b = 1; 
            c = 1;
            d = 1;
        elseif a == 2
            b = 2;
            c = 1;
            d = 1;
        elseif a == 3
            b = 2;
            c = 2;
            d = 1;
        elseif a == 4
            b = 1;
            c = 2;
            d = 1;
        elseif a == 5
            b = 1;
            c = 1;
            d = 2;
        elseif a == 6
            b = 2;
            c = 1;
            d = 2;
        elseif a == 7
            b = 2;
            c = 2;
            d = 2;
        elseif a == 8
            b = 1;
            c = 2;
            d = 2;
        end
        l_b      = PolyShape(deg_xi, b, xi, 0);
        l_b_xi   = PolyShape(deg_xi, b, xi, 1);
        l_c      = PolyShape(deg_eta, c, eta, 0);
        l_c_eta  = PolyShape(deg_eta, c, eta, 1);
        l_d      = PolyShape(deg_zeta, d, zeta, 0);
        l_d_zeta = PolyShape(deg_zeta, d, zeta, 1);
        val_xi = l_b_xi * l_c * l_d;
        val_eta = l_b * l_c_eta * l_d;
        val_zeta = l_b * l_c * l_d_zeta;

    % triquadratic hexahedral 27-nodes
    case 2
        if a == 1
            b = 1;
            c = 1;
            d = 1;
        elseif a == 2
            b = 3;
            c = 1;
            d = 1;
        elseif a == 3
            b = 3;
            c = 3;
            d = 1;
        elseif a == 4
            b = 1;
            c = 3;
            d = 1;
        elseif a == 5
            b = 1;
            c = 1;
            d = 3;
        elseif a == 6
            b = 3;
            c = 1;
            d = 3;
        elseif a == 7
            b = 3;
            c = 3;
            d = 3;
        elseif a == 8
            b = 1;
            c = 3;
            d = 3;
        elseif a == 9
            b = 2;
            c = 1;
            d = 1;
        elseif a == 10
            b = 3;
            c = 2; 
            d = 1;
        elseif a == 11
            b = 2;
            c = 3;
            d = 1;
        elseif a == 12
            b = 1;
            c = 2;
            d = 1;
        elseif a == 13
            b = 2;
            c = 1;
            d = 3;
        elseif a == 14
            b = 3;
            c = 2;
            d = 3;
        elseif a == 15
            b = 2; 
            c = 3;
            d = 3;
        elseif a == 16
            b = 1;
            c = 2;
            d = 3;
        elseif a == 17
            b = 1;
            c = 1;
            d = 2;
        elseif a == 18
            b = 3;
            c = 1;
            d = 2;
        elseif a == 19
            b = 3;
            c = 3;
            d = 2;
        elseif a == 20
            b = 1;
            c = 3;
            d = 2;
        elseif a == 21
            b = 2;
            c = 2;
            d = 1;
        elseif a == 22
            b = 2;
            c = 2;
            d = 3;
        elseif a == 23
            b = 2;
            c = 1;
            d = 2;
        elseif a == 24
            b = 2;
            c = 3;
            d = 2;
        elseif a == 25
            b = 1;
            c = 2;
            d = 2;
        elseif a == 26
            b = 3;
            c = 2;
            d = 2;
        elseif a == 27
            b = 2;
            c = 2;
            d = 2;
        end
        l_b      = PolyShape(deg_xi, b, xi, 0);
        l_b_xi   = PolyShape(deg_xi, b, xi, 1);
        l_c      = PolyShape(deg_eta, c, eta, 0);
        l_c_eta  = PolyShape(deg_eta, c, eta, 1);
        l_d      = PolyShape(deg_zeta, d, zeta, 0);
        l_d_zeta = PolyShape(deg_zeta, d, zeta, 1);
        val_xi   = l_b_xi * l_c * l_d;
        val_eta  = l_b * l_c_eta * l_d;
        val_zeta = l_b * l_c * l_d_zeta;

  
    otherwise
        error('Error: degree has to be 1, or 2.');
end