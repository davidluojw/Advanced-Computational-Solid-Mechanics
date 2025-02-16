function [K] = AssemblyK(deg,n_eq,n_en,n_ee, n_int,xi,weight,IEN,LM,n_el,x_coor,y_coor, mid_pt, n_ed, d, E, A, I)

% Allocate an empty stiffness matrix

K = sparse(n_eq, n_eq);
  
% Assembly the siffness matrix and load vector
for ee = 1 : n_el

    % element node
    node= [ x_coor(IEN(1, ee)), x_coor(IEN(n_en, ee));
                y_coor(IEN(1, ee)), y_coor(IEN(n_en, ee)) ];
    % displacement 
    uh = zeros(2,2);
    for aa  = 1:n_en
        for ii = 1:n_ed
            pp = n_ed * (aa - 1) + ii;
            PP = LM(pp, ee);
            if PP > 0
                if pp == 1
                    uh(1,1) = d(PP);
                elseif pp == 2
                    uh(2,1) = d(PP);
                elseif pp == 4
                    uh(1,2) = d(PP);
                elseif pp == 5
                    uh(2,2) = d(PP);
                end
            end
        end
    end
    % length changes due to displacement
    hh = elemLength(node, uh);

    % Allocate zero element stiffness matrix
    k_e = zeros(n_ee, n_ee);     % dimension of element stiffness is n_ee x n_ee
    
    % there is only one variable field in each element
    x_ele = zeros(n_en, 1);

    for aa = 1 : n_en
        AA = IEN(aa, ee);
        if AA <= mid_pt
            x_ele(aa) = y_coor(AA);
        else
            x_ele(aa) = x_coor(AA);
        end
    end
    
    % loop over quadrature points
    for ll = 1 : n_int

        x_l    = 0;
        dx_dxi = 0;
        
        for aa = 1 : n_en
            x_l    = x_l  + x_ele(aa) * PolyShape(deg, aa, xi(ll), 0);
            dx_dxi   = dx_dxi + x_ele(aa) * PolyShape(deg, aa, xi(ll), 1);
        end

        dxi_dx = 1.0 / dx_dxi;
          
        for aa = 1 : n_en
            for ii = 1: n_ed
                pp = n_ed * (aa - 1) + ii;
                
                for bb = 1 : n_en
                    for jj = 1 : n_ed
                        qq = n_ed * (bb - 1) + jj;
                        if ii == 1 && jj == 1
                            k_e(pp,qq) = k_e(pp,qq) + E * A * weight(ll) * PolyShape(deg, aa, xi(ll), 1) * PolyShape(deg, bb, xi(ll), 1) * dxi_dx;
                        end
                        if ii ~= 1 && jj ~= 1
                            pp_ind = (n_ed - 1) * (aa - 1) + ii - 1;
                            qq_ind = (n_ed - 1) * (bb - 1) + jj - 1;
                            k_e(pp,qq) = k_e(pp,qq) + E * I * weight(ll) * HermiteShape(pp_ind, xi(ll), 2, hh)...
                            * HermiteShape(qq_ind, xi(ll), 2, hh) * dxi_dx * dxi_dx * dxi_dx;
                        end
                    end
                end
            end
        end
    end
    % end of the quadrature loop

    % because of the displacement, the element angle changes
    angle = elemAngle(node, uh);

    % Compute the rotation matrix
    T = RotationMatrix(angle);

    % Trandform the stiffness matrix in terms of the global coordinates
    k_e = T * k_e * T';
    
    % Now we need to put element k and f into global K and F
    for aa = 1 : n_en
        for ii = 1 : n_ed
            pp = n_ed * (aa - 1) + ii;
            PP = LM(pp,ee);
            if PP > 0
                for bb = 1 : n_en
                    for jj = 1 : n_ed
                        qq = n_ed * (bb - 1) + jj;
                        QQ = LM(qq,ee);
                        if QQ > 0
                            K(PP,QQ) = K(PP,QQ) + k_e(pp,qq);
                        end
                    end
                end
            end
        end
    end
    
end




end