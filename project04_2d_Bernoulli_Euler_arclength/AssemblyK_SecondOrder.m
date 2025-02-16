function [K, Elem] = AssemblyK_SecondOrder(Model, Geometry, Solver, Elem, Elem_n, disp, disp_n)

% Allocate an empty stiffness matrix

K = sparse(Model.n_eq, Model.n_eq);
  
% Assembly the siffness matrix and load vector
for ee = 1 : Model.n_el

    % Element length
    L_1  = Elem_n(ee).eleLength;            % Initial load step element length
    hh   = elemLength(Elem_n(ee), disp);    % current element length
    D_L  = hh - L_1;                        % element axial deformation 

    % Rigid rotation increment
    rbr  = elemRigidRot(Elem_n(ee), disp, disp_n);
    % Absolute angle of the element
    angle = Elem_n(ee).eleAngle + rbr;

    % Update the angle 
    Elem(ee).eleAngle = angle;

    % Compute the rotation matrix
    T = RotationMatrix(angle);

    % deformation rotation, w.r.t initial load step
    deltadisp = disp - disp_n;
    r1 = deltadisp(Elem_n(ee).Node1.dof(3)) - rbr;
    r2 = deltadisp(Elem_n(ee).Node2.dof(3)) - rbr;

    % local coordinate, deformation increment, w.r.t initial load step
    dl = [0; 0 ; r1; D_L; 0; r2];

    % local coordinate, increment of internal force
    D_fl = Elem_n(ee).eleElasticK * dl;

    % local coordinate, total internal force
    fl = Elem_n(ee).eleFint + D_fl;

    % Update internal force
    Elem(ee).eleFint = fl;

    % Currertn iterate internal forces
    P  = Elem(ee).eleFint(4);
    M1 = Elem(ee).eleFint(3);
    M2 = Elem(ee).eleFint(6);

    
    % local coordinate, displacement & angle increment with respect to the
    
    % there is only one variable field in each element
    % x_ele = [0, hh];
    x_ele = [0; hh];


    % Allocate zero element stiffness matrix
    ke_ele = zeros(Model.n_ee, Model.n_ee);     % dimension of element stiffness is n_ee x n_ee
    kg_ele = zeros(Model.n_ee, Model.n_ee);
    
    % loop over quadrature points
    for ll = 1 : Model.n_int

        x_l    = 0;
        dx_dxi = 0;
        
        for aa = 1 : Model.n_en
            x_l    = x_l  + x_ele(aa) * PolyShape(Model.deg, aa, Solver.xi(ll), 0);
            dx_dxi   = dx_dxi + x_ele(aa) * PolyShape(Model.deg, aa, Solver.xi(ll), 1);
        end

        dxi_dx = 1.0 / dx_dxi;
          
        for aa = 1 : Model.n_en
            for ii = 1: Model.n_ed
                pp = Model.n_ed * (aa - 1) + ii;
                
                for bb = 1 : Model.n_en
                    for jj = 1 : Model.n_ed
                        qq = Model.n_ed * (bb - 1) + jj;
                        if ii == 1 && jj == 1
                            ke_ele(pp,qq) = ke_ele(pp,qq) + Model.E * Model.A * Solver.weight(ll) * PolyShape(Model.deg, aa, Solver.xi(ll), 1) *...
                                                            PolyShape(Model.deg, bb, Solver.xi(ll), 1) * dxi_dx;
                            kg_ele(pp,qq) = kg_ele(pp,qq) + P * Solver.weight(ll) * PolyShape(Model.deg, aa, Solver.xi(ll), 1) * ...
                                                            PolyShape(Model.deg, bb, Solver.xi(ll), 1) * dxi_dx;
                        end
                        if ii ~= 1 && jj ~= 1
                            pp_ind = (Model.n_ed - 1) * (aa - 1) + ii - 1;
                            qq_ind = (Model.n_ed - 1) * (bb - 1) + jj - 1;
                            ke_ele(pp,qq) = ke_ele(pp,qq) + Model.E * Model.I * Solver.weight(ll) * HermiteShape(pp_ind, Solver.xi(ll), 2, hh)...
                            * HermiteShape(qq_ind, Solver.xi(ll), 2, hh) * dxi_dx * dxi_dx * dxi_dx;
                            % kg_ele(pp,qq) = ke_ele(pp,qq) + Model.E * Model.I * Solver.weight(ll) * HermiteShape(pp_ind, Solver.xi(ll), 2, hh)...
                            % * HermiteShape(qq_ind, Solver.xi(ll), 2, hh) * dxi_dx * dxi_dx * dxi_dx;
                        end
                    end
                end
            end
        end
    end
    % end of the quadrature loop

    % generate the tangent matrix = elasctic matrix + geometric matrix
    kt_ele = ke_ele + kg_ele;

    % Store elastic stiffness matrix to be used in computation of internal forces
    Elem(ee).eleElasticK = ke_ele;
    Elem(ee).eleTangentK = kt_ele;

    % Trandform the stiffness matrix in terms of the global coordinates
   
    kt_ele = T * kt_ele * T';
    
    % Now we need to put element k and f into global K and F
    for aa = 1 : Model.n_en
        for ii = 1 : Model.n_ed
            pp = Model.n_ed * (aa - 1) + ii;
            PP = Geometry.LM(pp,ee);
            if PP > 0
                for bb = 1 : Model.n_en
                    for jj = 1 : Model.n_ed
                        qq = Model.n_ed * (bb - 1) + jj;
                        QQ = Geometry.LM(qq,ee);
                        if QQ > 0
                            K(PP,QQ) = K(PP,QQ) + kt_ele(pp,qq);
                        end
                    end
                end
            end
        end
    end
    
end




end