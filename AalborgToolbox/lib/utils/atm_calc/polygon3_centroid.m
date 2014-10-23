function c = polygon3_centroid(pol)
   
    proj_xy = [1 0 0;
               0 1 0;
               0 0 0];

    proj_xz = [1 0 0;
               0 0 0;
               0 0 1];
           
    proj_yz = [0 0 0;
               0 1 0;
               0 0 1];
   
    done = [0 0 0];
    pol_xy = proj_xy*pol;
    [tx, ty, e] = polygon2_centroid(pol_xy(1:2, :));
    cx = 0;
    cy = 0;
    cz = 0;
    if e == 0
        cx = tx;
        cy = ty;
        done(1) = 1;
        done(2) = 1;
    end
    pol_xz = proj_xz*pol;
    [tx, tz, e] = polygon2_centroid([pol_xz(1, :); pol_xz(3, :)]);
    if e == 0
        cx = tx;
        cz = tz;
        done(1) = 1;
        done(3) = 1;
    end
    if done(2) == 0 || done(3) == 0
        pol_yz = proj_yz*pol;
        [ty, tz, e] = polygon2_centroid([pol_yz(2, :); pol_yz(3, :)]);
        if e == 0
            if ~done(2)
                cy = ty;
            end
            if ~done(3)
                cz = tz;
            end
            done(2) = 1;
            done(3) = 1;
        end
    end
    if done(1) == 0
        cx = pol(1,1);
    end
    if done(2) == 0
        cy = pol(2,1);
    end    
    if done(3) == 0
        cz = pol(3,1);
    end
    c = [cx cy cz]';
end