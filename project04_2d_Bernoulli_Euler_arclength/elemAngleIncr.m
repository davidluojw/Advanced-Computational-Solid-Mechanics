function ang_incr = elemAngleIncr(node, uh_n, uh)

    % initial coordinate of the current load step = coordinate of the
    % previous load step
    x1_n = node(1,1) + uh_n(1,1);
    y1_n = node(2,1) + uh_n(2,1);
    x2_n = node(1,2) + uh_n(1,2);
    y2_n = node(2,2) + uh_n(2,2);

    % 当前荷载步初始单元轴向向量
    % Element axial direction of previous load step / initial current load
    % step
    vx_n = [x2_n-x1_n, y2_n-y1_n, 0];

    %当前迭代步更新后的新的节点坐标
    % updated coordinate of the current load step
    x1 = node(1,1) + uh(1,1);
    y1 = node(2,1) + uh(2,1);
    x2 = node(1,2) + uh(1,2);
    y2 = node(2,2) + uh(2,2);

    % 当前迭代步更新后的新的轴向向量
    % Elment axial direction of the current load step
    vx = [x2-x1, y2-y1, 0];

    % 角度增量
    % increment of the angle
    dir = cross(vx_n,vx);
    ang_incr = atan2(norm(dir),dot(vx_n,vx));

    % 角度转动方向
    % direction of the rotation
    s = sign(dir(3));
    ang_incr = s * ang_incr;
end