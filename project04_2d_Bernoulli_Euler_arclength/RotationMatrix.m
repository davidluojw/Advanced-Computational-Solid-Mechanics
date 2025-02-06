function T = RotationMatrix(node1, node2)
    % 计算单元方向向量
    dx = node2(1) - node1(1);
    dy = node2(2) - node1(2);
    L = sqrt(dx^2 + dy^2);
    
    % 计算方向余弦
    cos_theta = dx / L;
    sin_theta = dy / L;
    
    % 构造旋转矩阵
    T = [cos_theta,  sin_theta, 0, 0, 0, 0;
         -sin_theta, cos_theta, 0, 0, 0, 0;
         0,          0,         1, 0, 0, 0;
         0,          0,         0, cos_theta, sin_theta, 0;
         0,          0,         0, -sin_theta, cos_theta, 0;
         0,          0,         0, 0, 0, 1];
    
end
