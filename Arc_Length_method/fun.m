% out = [y,dy/dx]
function [y,dy]=fun(x)

    st=sin(pi/3);
    temp=1-2*st*x+x^2;
    y=(1/sqrt(temp)-1)*(st-x);
    dy=(st-x)^2/temp^(3/2)-(1/sqrt(temp)-1);
    
end