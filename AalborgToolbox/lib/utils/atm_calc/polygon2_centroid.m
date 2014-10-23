function [cx, cy, e] = polygon2_centroid(pol)

    A_v = [pol(1, end) pol(1, 1:end-1) pol(1, 1:end-1) pol(1, end)].*...
          [pol(2, 1:end-1) pol(2, end) -pol(2, end) -pol(2, 1:end-1)];
    A = sum(A_v)/2;
    
    if abs(A) < 1e-6
        e = 1;
        cx = 0;
        cy = 0;
    else
        e = 0;

        x_v = [pol(1, end) pol(1, 1:end-1)]+[pol(1, 1:end-1) pol(1, end)];
        cx = sum([x_v x_v].*A_v)/(6*A);

        y_v = [pol(2, end) pol(2, 1:end-1)]+[pol(2, 1:end-1) pol(2, end)];
        cy = sum([y_v y_v].*A_v)/(6*A);
    end
end