function ret = polygon_plane_division(pol, plane_n, plane_p)

    pol1 = [];
    pol2 = [];
    
    pols = [pol; pol(:,2:end) pol(:,1)];

    side = 0;
    for p = pols(:,:)
        [P, check] = plane_line_intersect(plane_n, plane_p, p(1:3,:), p(4:6,:));
        if check == 1
            pol1 = [pol1 P];
            pol2 = [pol2 P];
            side = ~side;
        end
        if side == 0
            pol1 = [pol1 p(4:6,:)];
        else
            pol2 = [pol2 p(4:6,:)];
        end
    end
    
    ret = {pol1, pol2};
    
end