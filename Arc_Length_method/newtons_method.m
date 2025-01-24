% 牛顿法
% input: 初始点x0, 载荷增量dy, 误差限tol
%output: 结束点[x1,y1]
function [x1,y1]=newtons_method(x0,dy,tol)
    
    xi=x0;
    [yi,dyi]=fun(xi);%初始点
    %dy=0.01;%载荷增量
    y1=yi+dy;%结束点

    Nmax=1000;num=1;
    while 1
        xj=xi+(y1-yi)/dyi;
        [yj,dyj]=fun(xj);
        xi=xj;yi=yj;dyi=dyj;
        num=num+1;
        if (abs(yi-y1)<tol) || (num==Nmax)
            break;
        end
    end
    x1=xi;y1=yi;

end