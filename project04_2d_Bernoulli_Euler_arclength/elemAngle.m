function angle = elemAngle(x_coor, y_coor, disp, AA, PP)
    % node1 and node2
    node1_old = [x_coor(AA(1)), y_coor(AA(1))];
    node2_old = [x_coor(AA(2)), y_coor(AA(2))];

    % extract current displacement from disp
    u1 = [disp(PP(1)), disp(PP(2))];
    u2 = [disp(PP(4)), disp(PP(5))];
    node1_new = node1_old + u1;
    node2_new = node2_old + u2;

    % axial direction
    axial_direction = node2_new - node1_new;

    % compute the new element length
    angle = atan2( axial_direction(2), axial_direction(1) );
end