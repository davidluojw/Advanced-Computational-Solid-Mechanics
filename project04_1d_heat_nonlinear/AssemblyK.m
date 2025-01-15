function [K] = AssemblyK(pp,n_eq,n_en,nqp,qp,wq,IEN,ID,nElem,uh,x_coor,fun_kappa,fun_dkappa)

% Allocate an empty stiffness matrix

K = sparse(n_eq, n_eq);

  
% Assembly the siffness matrix and load vector
for ee = 1 : nElem
    % Allocate zero element stiffness matrix
    k_ele = zeros(n_en, n_en);
    
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
            
            for bb = 1 : n_en
            Nb    = PolyBasis(pp, bb, 0, qp(qua));
            Nb_xi = PolyBasis(pp, bb, 1, qp(qua));
            k_ele(aa,bb) = k_ele(aa,bb) + wq(qua) * Na_xi * kappa * Nb_xi * dxi_dx;
            k_ele(aa,bb) = k_ele(aa,bb) + wq(qua) * Na_xi * dkappa * Nb * u_xi * dxi_dx;
            end
        end
    end
    % end of the quadrature loop
    
    % distribute the entries to the global stiffness matrix and global load vector
    for aa = 1 : n_en
        LM_a = ID( IEN(aa, ee) );
        if LM_a > 0
            for bb = 1 : n_en
                LM_b = ID( IEN(bb, ee) );
                if LM_b > 0
                  K(LM_a, LM_b) = K(LM_a, LM_b) + k_ele(aa, bb);
                end
            end
        end
    end
    
end

end