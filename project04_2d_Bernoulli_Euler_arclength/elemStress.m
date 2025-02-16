function sigma = elemStress(x_coor, y_coor, disp, deltadisp, AA, PP, E, Rigidrot)
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
    dv_dx = [0, HermiteShape(1,-1,1,L), HermiteShape(2,-1,1,L), ...
             0, HermiteShape(3,1,1,L), HermiteShape(4,1,1,L)] * d_e;

    sigma = E * [ du_dx + 0.5 * du_dx * du_dx + 0.5 * dv_dx * dv_dx;  % consider the Green-Lagrangian strain
                  0];  % neglect the higher-order term
                
end