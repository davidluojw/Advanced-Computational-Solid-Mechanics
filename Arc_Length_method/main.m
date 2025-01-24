%version - R2018b
%%
clc;clear;close all;

pt=12;%字体大小
tol=1e-10;%误差限
num=35;%点的个数
time=3;%绘图时间
dy=0.8/num;%载荷增量
x_newton=zeros(num,1);y_newton=x_newton;
dl=0.1;%弧长
x_arc=zeros(num,1);y_arc=x_arc;


x=0:0.01:3;
[y,~]=arrayfun(@fun,x);
h=figure;figset([20,7,0.3,0.3]);
subplot(1,2,1);plot(x,y,'k','linewidth',1.5);hold on;grid on;
set(gca,'ticklabelinterpreter','latex','fontsize',pt);
xlabel('$a$','interpreter','latex','fontsize',pt);
ylabel('$\lambda$','interpreter','latex','fontsize',pt,'Rotation',0);
title('Newton''s Method','interpreter','latex','fontsize',pt);

subplot(1,2,2);plot(x,y,'k','linewidth',1.5);hold on;grid on;
set(gca,'ticklabelinterpreter','latex','fontsize',pt);
xlabel('$a$','interpreter','latex','fontsize',pt);
ylabel('$\lambda$','interpreter','latex','fontsize',pt,'Rotation',0);
title('Arc Length Method','interpreter','latex','fontsize',pt);

for i=2:num
    [x_newton(i),y_newton(i)]=newtons_method(x_newton(i-1),dy,tol);
    subplot(1,2,1);
    plot(x_newton(i),y_newton(i),'ob','linewidth',1.5,'markersize',7);
    
    [x_arc(i),y_arc(i)]=arc_method(x_arc(i-1),dl,tol);
    subplot(1,2,2);
    plot(x_arc(i),y_arc(i),'ob','linewidth',1.5,'markersize',7);
    pause(time/num);
end

%%
%设置图片尺寸与位置
function figset(parameter1,parameter2)

    %电脑屏幕尺寸
    set(0,'units','centimeters')
    cm_ss=get(0,'screensize');
    W=cm_ss(3);%电脑屏幕长度，单位cm
    H=cm_ss(4);%电脑屏幕宽度，单位cm

    %设置figure在screen中的比例位置与大小
    temp1=[parameter1(3),parameter1(4),parameter1(1)/W,parameter1(2)/H];
    set(gcf,'units','normalized','position',temp1);
    if nargin==2
        %设置axis在figure中的比例位置与大小
        temp2=[parameter2(3),parameter2(4),parameter2(1),parameter2(2)];
        set(gca,'position',temp2);
    end

end