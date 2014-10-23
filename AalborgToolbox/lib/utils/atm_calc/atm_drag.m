function out = atm_drag(input)

    rho = input(1);
    Cd = input(2);
    com = input(3:5);
    v = input(6:8);
    r = input(9:11);
    q1 = input(12);
    q2 = input(13);
    q3 = input(14);
    q4 = input(15);
    
    w = [0 0 0.7292e-4]';
    cv = -cross(w, r);
    
    v = v+cv;
    
    v_u = v/norm(v);

    R_ib = [q1^2-q2^2-q3^2+q4^2, 2*(q1*q2+q3*q4), 2*(q1*q3-q2*q4);
            2*(q1*q2-q3*q4), -q1^2+q2^2-q3^2+q4^2, 2*(q2*q3+q1*q4);
            2*(q1*q3+q2*q4), 2*(q2*q3-q1*q4), -q1^2-q2^2+q3^2+q4^2];
                
    small_plate = [0.0 -0.05 -0.05;
                   0.0  0.05 -0.05;
                   0.0  0.05  0.05;
                   0.0 -0.05  0.05]';

    large_plate = [0.1035 0.0 -0.05;
                  -0.1035 0.0 -0.05;
                  -0.1035 0.0  0.05;
                   0.1035 0.0  0.05]';

    antenna = [0.0    0.005 0.0;
               0.0   -0.005 0.0;
               -0.15 -0.005 0.0;
               -0.15  0.005 0.0]';

    plate_front    = translate(rotate2(small_plate, [  0   0   0]), [ 0.1035 0 0]);
    plate_back     = translate(rotate2(small_plate, [  0 180   0]), [-0.1035 0 0]);
    plate_left     = translate(rotate2(large_plate, [  0   0 180]), [0 -0.05 0]);
    plate_right    = translate(rotate2(large_plate, [  0   0   0]), [0  0.05 0]);
    plate_top      = translate(rotate2(large_plate, [ 90   0   0]), [0  0 -0.05]);
    plate_bottom   = translate(rotate2(large_plate, [-90   0   0]), [0  0  0.05]);
    antenna_top_f    = translate(rotate2(antenna, [  0 -45 0]),  [-0.1035  0     0.05]);
    antenna_bottom_f = translate(rotate2(antenna, [180  45 0]),  [-0.1035  0    -0.05]);
    antenna_left_f   = translate(rotate2(antenna, [-90  0 -45]), [-0.1035 -0.05  0]);
    antenna_right_f  = translate(rotate2(antenna, [ 90  0  45]), [-0.1035  0.05  0]);
    antenna_top_b    = translate(rotate2(antenna, [180 -45 0]),  [-0.1035  0     0.05]);
    antenna_bottom_b = translate(rotate2(antenna, [  0  45 0]),  [-0.1035  0    -0.05]);
    antenna_left_b   = translate(rotate2(antenna, [ 90  0 -45]), [-0.1035 -0.05  0]);
    antenna_right_b  = translate(rotate2(antenna, [-90  0  45]), [-0.1035  0.05  0]);


    satellite_b = {plate_front, plate_back, plate_left, plate_right, ...
                    plate_top, plate_bottom, antenna_top_f, antenna_bottom_f,...
                    antenna_left_f, antenna_right_f, antenna_top_b,...
                    antenna_bottom_b, antenna_left_b, antenna_right_b};
                

    mass_center = com;

    for i = 1:length(satellite_b)
        satellite_b{i} = polygon_backface_check(satellite_b{i}, -R_ib*v_u);
    end
    satellite_b = flatten(satellite_b);
    
    F = 0;
    T = 0;
    for i = 1:length(satellite_b)
        c = polygon3_centroid(satellite_b{i})-mass_center;
        A = polygon3_area(satellite_b{i}, R_ib*v_u);
        Ft = -v_u*(rho*Cd*A*norm(v)^2)/2;
        Tt = cross(c, R_ib*Ft);

        F = F + Ft;
        T = T + Tt;
    end
    
    out = [T; F];
    
end
