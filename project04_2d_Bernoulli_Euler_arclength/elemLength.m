function L = elemLength(Elem, disp)
    % Element nodes
    Node1 = Elem.Node1;
    Node2 = Elem.Node2;

    % Nodal displacements
    dx1 = disp(Node1.dof(1));
    dy1 = disp(Node1.dof(2));
    dx2 = disp(Node2.dof(1));
    dy2 = disp(Node2.dof(2));

    % New nodal coordinates
    x1 = Node1.x + dx1;
    y1 = Node1.y + dy1;
    x2 = Node2.x + dx2;
    y2 = Node2.y + dy2;

    % compute the new element length
    L = norm( [x2 - x1,  y2 - y1] );
end