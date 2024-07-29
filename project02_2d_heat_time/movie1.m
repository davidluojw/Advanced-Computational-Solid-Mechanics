% 定义网格大小
nx = 50;
ny = 50;
nz = 20;

% 生成网格点
[x, y, z] = meshgrid(linspace(-5, 5, nx), linspace(-5, 5, ny), linspace(-2, 2, nz));

% 生成模拟温度数据 (例如，高斯分布)
temperature = exp(-0.1 * (x.^2 + y.^2 + z.^2));

% 创建视频写入对象
v = VideoWriter('temperature_field_animation.avi');
v.FrameRate = 10; % 设置帧率
open(v);

% 创建图形窗口
figure;

% 定义切片位置
slice_z = linspace(-2, 2, 10);

% 初始化帧数组
frames(length(slice_z)) = struct('cdata', [], 'colormap', []);

% 循环生成每一帧
for i = 1:length(slice_z)
    % 清空当前图形
    clf;
    
    % 创建切片图
    slice(x, y, z, temperature, [], [], slice_z(i));
    shading interp;
    colorbar;
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title(['3D Temperature Field Slice at z = ', num2str(slice_z(i))]);
    
    % 获取当前帧
    frame = getframe(gcf);
    
    % 保存帧到数组
    frames(i) = frame;
    
    % 写入帧到视频
    writeVideo(v, frame);
end

% 关闭视频写入对象
close(v);

% 播放动画
figure;
movie(frames, 1, v.FrameRate);
