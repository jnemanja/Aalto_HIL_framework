function n = polygon_normal(pol)

    n = sum(cross(pol, [pol(:,2:end) pol(:,1)]),2);
    
end