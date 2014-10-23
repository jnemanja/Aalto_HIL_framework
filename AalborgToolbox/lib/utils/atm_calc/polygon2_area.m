function A = polygon2_area(pol)

    A_v = [pol(1, end) pol(1, 1:end-1) pol(1, 1:end-1) pol(1, end)].*...
          [pol(2, 1:end-1) pol(2, end) -pol(2, end) -pol(2, 1:end-1)];
    A = abs(sum(A_v)/2);

end