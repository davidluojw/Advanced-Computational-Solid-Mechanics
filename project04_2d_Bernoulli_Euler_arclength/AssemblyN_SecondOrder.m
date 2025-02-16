function [N, Elem] = AssemblyN_SecondOrder(Model, Geometry, Solver, Elem, Elem_n, disp, disp_n)
% Allocate an empty load vector
N = zeros(Model.n_eq, 1);
  
% Assembly of K and F
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


    % 内力从局部到全局的转换
    n_e = T * fl;


    % Now we need to put element internal force n_e into global N
    for aa = 1 : Model.n_en
        for ii = 1 : Model.n_ed
            pp = Model.n_ed * (aa - 1) + ii;
            PP = Geometry.LM(pp,ee);
            if PP > 0
                N(PP) = N(PP) + n_e(pp);
            end
        end
    end

    
end

end