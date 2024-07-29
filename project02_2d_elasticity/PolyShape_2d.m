function val = PolyShape_2d(deg_xi, deg_eta, a, xi, eta)

if deg_xi == deg_eta
    degree = deg_xi;
else
    error('Error: the two degree of xi and eta must be identical.');
end

switch degree
    % bilinear quadraliteral 4-nodes
    case 1
        if a == 1
            b = 1; 
            c = 1;
        elseif a == 2
            b = 2;
            c = 1;
        elseif a == 3
            b = 2;
            c = 2;
        elseif a == 4
            b = 1;
            c = 2;
        end
        l_b = PolyShape(deg_xi, b, xi, 0);
        l_c = PolyShape(deg_eta, c, eta, 0);
        val = l_b * l_c;

    % biquadratic quadraliteral 9-nodes
    case 2
        if a == 1
            b = 1;
            c = 1;
        elseif a == 2
            b = 3;
            c = 1;
        elseif a == 3
            b = 3;
            c = 3;
        elseif a == 4
            b = 1;
            c = 3;
        elseif a == 5
            b = 2;
            c = 1;
        elseif a == 6
            b = 3;
            c = 2;
        elseif a == 7
            b = 2;
            c = 3;
        elseif a == 8
            b = 1;
            c = 2;
        elseif a == 9
            b = 2;
            c = 2;
        end
        l_b = PolyShape(deg_xi, b, xi, 0);
        l_c = PolyShape(deg_eta, c, eta, 0);
        val = l_b * l_c;

    % bicubic quadraliteral 16-nodes
    case 3
        if a == 1
            b = 1;
            c = 1;
        elseif a == 2
            b = 4;
            c = 1;
        elseif a == 3
            b = 4;
            c = 4;
        elseif a == 4
            b = 1;
            c = 4;
        elseif a == 5
            b = 2;
            c = 1;
        elseif a == 6
            b = 4;
            c = 2;
        elseif a == 7
            b = 3;
            c = 4;
        elseif a == 8
            b = 1;
            c = 3;
        elseif a == 9
            b = 3;
            c = 1;
        elseif a == 10
            b = 4;
            c = 3;
        elseif a == 11
            b = 2;
            c = 4;
        elseif a == 12
            b = 1;
            c = 2;
        elseif a == 13
            b = 2;
            c = 2;
        elseif a == 14
            b = 3;
            c = 2;
        elseif a == 15
            b = 3;
            c = 3;
        elseif a == 16
            b = 2;
            c = 3;
        end
        l_b = PolyShape(deg_xi, b, xi, 0);
        l_c = PolyShape(deg_eta, c, eta, 0);
        val = l_b * l_c;

    otherwise
        error('Error: degree has to be 1, 2, or 3.');
end