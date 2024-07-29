% 创建一个新的图形窗口
figure;

% 定义 x 轴数据
x = linspace(0, 2*pi, 100);

% 定义动画的帧数
numFrames = 100;

% 初始化存储帧的结构体数组
F(numFrames) = struct('cdata', [], 'colormap', []);

% 设置轴范围
axis([0 2*pi -1 1]);

% 标题
title('Sine Wave Animation');
xlabel('x');
ylabel('sin(x)');

for k = 1:numFrames
    % 更新 y 轴数据，制造动画效果
    y = sin(x - 2*pi*k/numFrames);
    
    % 绘制数据
    plot(x, y, 'b', 'LineWidth', 2);
    
    % 捕获当前图形窗口的内容为帧
    F(k) = getframe(gcf);
    
    % 暂停一会儿以控制帧速率
    pause(0.05);
end

% 播放动画
movie(gcf, F, 2, 15); % 播放2次，15帧每秒

% 创建视频写入对象
v = VideoWriter('sine_wave_animation.avi');
open(v);

% 写入每一帧到视频
for k = 1:numFrames
    writeVideo(v, F(k));
end

% 关闭视频文件
close(v);


%-------------------------------------------------------------------------
% 制作三维动画
% 创建一个新的图形窗口
figure;

% 定义 x 和 y 轴数据
[x, y] = meshgrid(linspace(-2, 2, 50), linspace(-2, 2, 50));

% 定义动画的帧数
numFrames = 100;

% 初始化存储帧的结构体数组
F(numFrames) = struct('cdata', [], 'colormap', []);

% 设置轴范围
axis([-2 2 -2 2 -2 2]);

% 标题
title('3D Surface Animation');
xlabel('x');
ylabel('y');
zlabel('z');

for k = 1:numFrames
    % 更新 z 轴数据，制造动画效果
    z = sin(sqrt(x.^2 + y.^2) - 2*pi*k/numFrames);
    
    % 绘制数据
    surf(x, y, z);
    
    % 捕获当前图形窗口的内容为帧
    F(k) = getframe(gcf);
    
    % 暂停一会儿以控制帧速率
    pause(0.05);
end
% 播放动画
movie(gcf, F, 2, 15); % 播放2次，15帧每秒
% 创建视频写入对象
v = VideoWriter('3d_surface_animation.avi');
open(v);

% 写入每一帧到视频
for k = 1:numFrames
    writeVideo(v, F(k));
end

% 关闭视频文件
close(v);



