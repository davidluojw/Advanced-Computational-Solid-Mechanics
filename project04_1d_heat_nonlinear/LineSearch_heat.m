function [s] = LineSearch_heat(G0,G,Deltad,d,stol,omega_r,g,pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h)
s = 1.0;
% Line Search
% find a s such that abs(G) <= stol * abs(G0)
% if abs(G) > stol * abs(G0)
linmax = 10;
smax = 16;
sb = 0;
sa = s;
Gb = G0;
Ga = G;

% temporary variables
d_temp = d;

% find bracket on zero
% if Ga * Gb > 0, [sb,sa] no zero point
% then enlarge the bracket to [sa, 2sa]
% compute the corresponding R_temp and Gb, Ga 
% in new bracket [sa,2sa]
while (Ga * Gb > 0 && sa < smax)
    sb = sa;
    sa = 2*sa;
    Gb = Ga;
    d_temp = d_temp + sa * Deltad;
    uh = [ d_temp ; g(omega_r) ];
    F_temp = AssemblyF(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);
    Ga = Deltad' * F_temp;
end

step = sa; 
G = Ga;
lin_step = 0;
% now we have already found the bracket [sb,sa]
% which contains the zero point. Ga * Gb < 0

% Illinois algorithm to find zero
% while still abs(G) > stol * abs(G0), criteria is not satisfied
while (lin_step <= linmax && Ga * Gb < 0 && ...
        abs(G) > stol * abs(G0) ) %|| abs(sb-sa) > stol * 0.5 * (sb+sa)))
    step = sa - Ga * (sa - sb) / (Ga - Gb);   % linear interpolation
    d_temp = d_temp + step * Deltad;
    uh = [ d_temp ; g(omega_r) ];
    F_temp = AssemblyF(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h);
    G = Deltad' * F_temp;

    if G * Ga > 0
        Gb = 0.5 * Gb;
    else
        sb = sa;
        Gb = Ga;
    end
    sa = step;
    Ga = G;
    lin_step = lin_step + 1;
end
s = step;
% update the new updates
% d = d - Deltad;
% d = d + s * Deltad;
    
% end
end