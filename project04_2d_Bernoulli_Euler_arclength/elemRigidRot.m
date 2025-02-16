function ang_incr = elemRigidRot(Elem, disp, disp_n)
    % Element nodes
    Node1 = Elem.Node1;
    Node2 = Elem.Node2;

    % Initial load step Nodal displacements
    dx1_0 = disp_n(Node1.dof(1));
    dy1_0 = disp_n(Node1.dof(2));
    dx2_0 = disp_n(Node2.dof(1));
    dy2_0 = disp_n(Node2.dof(2));

    % Initial load step Nodal coordinates
    x1_0 = Node1.x + dx1_0;
    y1_0 = Node1.y + dy1_0;
    x2_0 = Node2.x + dx2_0;
    y2_0 = Node2.y + dy2_0;

    % Initial load step element axial direction
    vx_0 = [x2_0 - x1_0, y2_0 - y1_0, 0];

    % Current iterate step Nodal displacements
    dx1_1 = disp(Node1.dof(1));
    dy1_1 = disp(Node1.dof(2));
    dx2_1 = disp(Node2.dof(1));
    dy2_1 = disp(Node2.dof(2));

    % Current iterate step Nodal coordinates
    x1_1 = Node1.x + dx1_1;
    y1_1 = Node1.y + dy1_1;
    x2_1 = Node2.x + dx2_1;
    y2_1 = Node2.y + dy2_1;

    % Current iterate step element axial direction
    vx_1 = [x2_1 - x1_1, y2_1 - y1_1, 0];

    % Increment of the angle
    dir = cross(vx_0,vx_1);
    ang_incr = atan2(norm(dir),dot(vx_0,vx_1));

    % direction the rotation
    s = sign(dir(3));
    ang_incr = s * ang_incr;
end