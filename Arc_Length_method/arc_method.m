% 弧长法
% input: 初始点x0, 弧长dl, 误差限tol
%output: 结束点[x1,y1]
function [x1,y1]=arc_method(x0,dl,tol)

    %tol=1e-10;%误差限
    [y0,dy0]=fun(x0);%初始点
    xi=x0;yi=y0;dyi=dy0;
    %dl=0.01;%弧长

    Nmax=1000;num=1;
    while 1
        temp_t=yi-dyi*xi-y0;
        temp_a=1+dyi^2;
        temp_b=2*dyi*temp_t-2*x0;
        temp_c=temp_t^2+x0^2-dl^2;
        delta=temp_b^2-4*temp_a*temp_c;
        if delta<0
            disp('delta<0');
            break;
        end
        xj1=(-temp_b+sqrt(delta))/2/temp_a;
        xj2=(-temp_b-sqrt(delta))/2/temp_a;
        xj=max(xj1,xj2);%取增大项
        [yj,dyj]=fun(xj);
        num=num+1;
        if (abs(yj-yi)<tol) || (num==Nmax)
            %disp(['err=',num2str(abs(yj-yi))]);
            %disp(['radius=',num2str(sqrt((xj-x0)^2+(yj-y0)^2))]);
            break;
        end
        xi=xj;yi=yj;dyi=dyj;
    end
    x1=xi;y1=yi;
    
end