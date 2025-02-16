function Fint = elemFint(x_coor, y_coor, disp, deltadisp, AA, PP, E, I, A, Rigidrot)
    % element dofs
    node1_old = [x_coor(AA(1)), y_coor(AA(1))];
    node2_old = [x_coor(AA(2)), y_coor(AA(2))];

    % extract current displacement from disp
    u1 = [disp(PP(1)), disp(PP(2))];
    u2 = [disp(PP(3)), disp(PP(4))];
    node1_new = node1_old + u1;
    node2_new = node2_old + u2;

    % compute the new element length
    L = norm( node2_new - node1_new );

    % element dofs
    d_e = deltadisp(PP(1:6));
    d_e(3) = d_e(3) - Rigidrot;
    d_e(6) = d_e(6) - Rigidrot;
    du_dx = [PolyShape(1,1,-1,0), 0, 0, PolyShape(1,2,1,0), 0, 0] * d_e;
    dv_dxdx = [0, HermiteShape(1,-1,2,L), HermiteShape(2,-1,2,L), ...
               0, HermiteShape(3,1,2,L), HermiteShape(4,1,2,L)] * d_e;
    dv_dxdxdx = [0, HermiteShape(1,-1,3,L), HermiteShape(2,-1,3,L), ...
                 0, HermiteShape(3,1,3,L), HermiteShape(4,1,3,L)] * d_e;
    P = E * A * du_dx;
    Q = E * I * dv_dxdxdx;
    M = E * I * dv_dxdx;
    Fint = [P, Q, M]';

end