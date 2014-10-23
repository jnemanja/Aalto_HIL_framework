function res = polygon_backface_check(pol, v)

    if dot(polygon_normal(pol), v)<0 
        res = pol;
    else
        res = zeros(3, 4);
    end

end