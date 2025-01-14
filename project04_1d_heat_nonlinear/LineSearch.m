function [s] = LineSearch(G0,G,d,R1,R2,Deltad,n,stol)
s = 1;
% Line Search
% find a s such that abs(G) <= stol * abs(G0)
if abs(G) > stol * abs(G0)
    linmax = 10;
    smax = 16;
    sb = 0;
    sa = s;
    Gb = G0;
    Ga = G;

    % temporary variables
    d_temp = d;
    R_temp = zeros(2,1);

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
        R_temp(1) = R1(d_temp(1),d_temp(2),n);
        R_temp(2) = R2(d_temp(1),d_temp(2),n);
        Ga = Deltad' * R_temp;
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
        R_temp(1) = R1(d_temp(1),d_temp(2),n);
        R_temp(2) = R2(d_temp(1),d_temp(2),n);
        G = Deltad' * R_temp;

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
    
end
end