% 定义函数及其导数
f = @(x) x^3 - x - 1;
df = @(x) 3*x^2 - 1;
 
% 初始近似解与精度要求
x0 = 1.5;
E_target = 1e-6;
 
% 迭代过程
x = x0;
E = Inf;
iter = 0;
while E >= E_target
    x_prev = x;
    x = x - f(x)/df(x);
    E = abs((x - x_prev)/x);
    iter = iter + 1;
end
 
% 输出结果
disp(['近似解：', num2str(x)]);
disp(['迭代次数：', num2str(iter)]);
disp(['收敛速度：', num2str(E)]);
