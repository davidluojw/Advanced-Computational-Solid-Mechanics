function [F] = AssemblyF(Model, Geometry, BC, Solver)
% Allocate an empty load vector
F = zeros(Model.n_eq, 1);
  
% Assembly of K and F
for ee = 1 : Model.n_el

    % element node
    node1 = [ Geometry.x_coor(Geometry.IEN(1, ee)) ,   Geometry.y_coor(Geometry.IEN(1, ee))    ];
    node2 = [ Geometry.x_coor(Geometry.IEN(Model.n_en, ee)), Geometry.y_coor(Geometry.IEN(Model.n_en, ee)) ];
    % element length
    dx = node2(1) - node1(1);
    dy = node2(2) - node1(2);
    hh = sqrt(dx^2 + dy^2);


    f_e = zeros(Model.n_ee, 1);            % element force nodes

    % there is only one variable field in each element
    x_ele = zeros(Model.n_en, 1);

    for aa = 1:Model.n_en
        AA = Geometry.IEN(aa, ee);
        if AA <= Geometry.mid_pt
            x_ele(aa) = Geometry.y_coor(AA);
        else
            x_ele(aa) = Geometry.x_coor(AA);
        end
    end

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
                if ii == 1
                    f_e(pp) = f_e(pp) + Solver.weight(ll) * PolyShape(Model.deg, aa, Solver.xi(ll), 0) * BC.f(x_l) * dx_dxi;
                else
                    pp_ind = (Model.n_ed - 1) * (aa - 1) + ii - 1;
                    f_e(pp) = f_e(pp) + Solver.weight(ll) * HermiteShape(pp_ind, Solver.xi(ll), 0, hh) * BC.f(x_l) * dx_dxi;
                end
                
            end
        end
    end



    % Now we need to put element k and f into global K and F
    for aa = 1 : Model.n_en
        for ii = 1 : Model.n_ed
            pp = Model.n_ed * (aa - 1) + ii;
            PP = Geometry.LM(pp,ee);
            if PP > 0
                F(PP) = F(PP) + f_e(pp);
            end
        end
    end

    % if node2(1) == 24
    %     F(Geometry.LM(2, Geometry.IEN(2, ee))) = BC.F_applied;
    % end


    for aa = 1: Model.n_en
        AA = Geometry.IEN(aa, ee);
        if AA == floor(Model.n_np / 2) + 2
            F(Geometry.LM(2, AA)) = BC.F_applied;
        end
    end
end

end