function [N] = AssemblyN(deg,n_eq,n_en,n_ee, n_int,xi,weight,IEN,LM,n_el,x_coor,y_coor, mid_pt, n_ed, d, deltad, E, A, I, elemLength_UL, elemAngle_UL, elemFint_UL)
% Allocate an empty load vector
N = zeros(n_eq, 1);
  
% Assembly of K and F
for ee = 1 : n_el

    % current iterate step element length
    hh = elemLength_UL(ee, 2);      % updated element length

    % length changes due to displacement
    % element node coordinates, based on the previous load step
    node= [ x_coor(IEN(1, ee)), x_coor(IEN(n_en, ee));
            y_coor(IEN(1, ee)), y_coor(IEN(n_en, ee)) ];
    % displacement, based on the current iterate step
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
    % because of the displacement, the element length changes
    L_c = elemLength(node, uh);

    % total changes of length, from previous load step to current iterate
    % step
    deltaL = L_c - hh;

    % increment of the rigid body rotation, with respect to the previous
    % load step
    d_n = d - deltad;
    % displacement of the previous load step
    uh_n = zeros(2,2);
    for aa  = 1:n_en
        for ii = 1:n_ed
            pp = n_ed * (aa - 1) + ii;
            PP = LM(pp, ee);
            if PP > 0
                if pp == 1
                    uh_n(1,1) = d_n(PP);
                elseif pp == 2
                    uh_n(2,1) = d_n(PP);
                elseif pp == 4
                    uh_n(1,2) = d_n(PP);
                elseif pp == 5
                    uh_n(2,2) = d_n(PP);
                end
            end
        end
    end
    % increment of the rigid body rotation
    rigidbody_rotation = elemAngleIncr(node, uh_n, uh);
    % updated angle, with respect to the previous load step
    angle = elemAngle_UL(ee, 2) + rigidbody_rotation;

    update_angle = 1;
    if (update_angle)
        elemAngle_UL(ee,1) = angle;
    end

    % deformation rotation, with respect to the previous load step
    % configuration
    angle_deform = zeros(2, 1);
    for aa = 1:n_en
        for ii = 1:n_ed
            pp = n_ed * (aa - 1) + ii;
            PP = LM(pp, ee);
            if PP > 0
                if pp == 3
                    angle_deform(1) = deltad(PP) - rigidbody_rotation;
                elseif pp == 6
                    angle_deform(2) = deltad(PP) - rigidbody_rotation;
                end
            end
        end
    end

    % local coordinate, displacement & angle increment with respect to the
    % previous load step
    d_ele = [0; 0; angle_deform(1); deltaL; 0; angle_deform(2)];
    
    % local coordinate, there is only one variable field in each element
    x_ele = zeros(n_en, 1);

    for aa = 1:n_en
        AA = IEN(aa, ee);
        if AA <= mid_pt
            x_ele(aa) = y_coor(AA);
        else
            x_ele(aa) = x_coor(AA);
        end
    end

    
    % local coordinate, element internal force increment
    deltan_e = zeros(n_ee, 1);

    % loop over quadrature points
    for ll = 1 : n_int

        x_l    = 0;
        dx_dxi = 0;
        % d_l    = zeros(n_ee, 1);
        du_dxi  = 0;
        dv_dxixi = 0;
        
        for aa = 1 : n_en
            x_l    = x_l    + x_ele(aa) * PolyShape(deg, aa, xi(ll), 0);
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(deg, aa, xi(ll), 1);
            
            for ii = 1:n_ed
                pp = n_ed * (aa - 1) + ii;
                if ii == 1
                    du_dxi = du_dxi + d_ele(pp) * PolyShape(deg, aa, xi(ll), 1);
                else
                    pp_ind = (n_ed - 1) * (aa - 1) + ii - 1;
                    dv_dxixi = dv_dxixi + d_ele(pp) * HermiteShape(pp_ind, xi(ll), 2, hh);
                end
            end
            
        end

        dxi_dx = 1.0 / dx_dxi;
    
        for aa = 1 : n_en
            for ii = 1: n_ed
                pp = n_ed * (aa - 1) + ii;
                if ii == 1
                    % f_e(pp) = f_e(pp) + lambda * weight(ll) * PolyShape(deg, aa, xi(ll), 0) * f(x_l) * dx_dxi;
                    deltan_e(pp) = deltan_e(pp) + E * A * weight(ll) * PolyShape(deg, aa, xi(ll), 1) * du_dxi * dxi_dx;
                else
                    pp_ind = (n_ed - 1) * (aa - 1) + ii - 1;
                    % f_e(pp) = f_e(pp) + weight(ll) * HermiteShape(pp_ind, xi(ll), 0, hh) * f(x_l) * dx_dxi;
                    deltan_e(pp) = deltan_e(pp) + E * I * weight(ll) * HermiteShape(pp_ind, xi(ll), 2, hh) * dv_dxixi * dxi_dx * dxi_dx * dxi_dx;
                end
                
            end
        end
    end

    % local coordinates, total element internal force
    % update element internal force
    elemFint_UL(ee, 1:n_ee) = elemFint_UL(ee, (n_ee+1):end) + deltan_e';

    % Compute the rotation matrix
    T = RotationMatrix(angle);

    % global coordinate, element internal force
    n_e = T * elemFint_UL(ee, 1:n_ee)';

    % Now we need to put element internal force n_e into global N
    for aa = 1 : n_en
        for ii = 1 : n_ed
            pp = n_ed * (aa - 1) + ii;
            PP = LM(pp,ee);
            if PP > 0
                N(PP) = N(PP) + n_e(pp);
                % for bb = 1 : n_en
                %     for jj = 1 : n_ed
                %         qq = n_ed * (bb - 1) + jj;
                %         QQ = LM(qq,ee);
                %         if QQ <= 0
                %             F(PP) = F(PP) - k_e(pp,qq) * g1 - k_e(pp,qq) * g2;
                %         end
                %     end
                % end
            end
        end
    end

    
end

end