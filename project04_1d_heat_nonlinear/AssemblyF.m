function [F] = AssemblyF(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa,f,h)
% Allocate an empty load vector
F = zeros(n_eq, 1);
  
% Assembly the siffness matrix and load vector
for ee = 1 : nElem
    % Allocate zero  element load vector
    f_ele = zeros(n_en, 1);
    
    x_ele = zeros(n_en, 1);
    d_ele = zeros(n_en, 1);
    for aa = 1 : n_en
        x_ele(aa) = x_coor( IEN(aa, ee) );
        u_ele(aa) = uh( IEN(aa, ee) );
    end
    
    for qua = 1 : nqp
        % geometrical mapping
        x_qua    = 0.0;
        dx_dxi   = 0.0;
        u_qua    = 0.0;
        u_xi     = 0.0;
        for aa = 1 : n_en
            x_qua    = x_qua  + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            dx_dxi   = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
            u_qua    = u_qua  + u_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            u_xi     = u_xi   + u_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
        end
        dxi_dx = 1.0 / dx_dxi;
          
        kappa = fun_kappa( u_qua );
        dkappa = fun_dkappa( u_qua );
          
        for aa = 1 : n_en
            Na    = PolyBasis(pp, aa, 0, qp(qua));
            Na_xi = PolyBasis(pp, aa, 1, qp(qua));
            f_ele(aa) = f_ele(aa) + wq(qua) * Na * f(x_qua) * dx_dxi;
            f_ele(aa) = f_ele(aa) - wq(qua) * Na_xi * kappa * u_xi * dxi_dx;
           
        end
    end
    % end of the quadrature loop
    
    % distribute the entries to the global stiffness matrix and global load vector
    for aa = 1 : n_en
        LM_a = ID( IEN(aa, ee) );
        if LM_a > 0
            F(LM_a) = F(LM_a) + f_ele(aa);
        end
    end
    
    % Modify the load vector by the Natural BC
    % Note: for multi-dimensional cases, one needs to perform line or
    % surface integration for the natural BC.
    if ee == 1
        F( ID(IEN(1, ee)) ) = F( ID(IEN(1, ee)) ) + h( x_coor(IEN(1,ee)));
    end
end
end