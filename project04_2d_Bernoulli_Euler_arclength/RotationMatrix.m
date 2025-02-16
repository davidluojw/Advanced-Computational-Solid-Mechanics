function T = RotationMatrix(angle)
    
    % 计算方向余弦
    cos_theta = cos(angle);
    sin_theta = sin(angle);
    
    % 构造旋转矩阵
    T = [cos_theta, -sin_theta, 0, 0, 0, 0;
         sin_theta, cos_theta, 0, 0, 0, 0;
         0,          0,         1, 0, 0, 0;
         0,          0,         0, cos_theta, -sin_theta, 0;
         0,          0,         0, sin_theta, cos_theta, 0;
         0,          0,         0, 0, 0, 1];
    
end
