function A = polygon3_area(pol, n)

    A_v = cross([pol(:, end) pol(:, 1:end-1)], [pol(:, 1:end-1) pol(:, end)]);
    A = abs(dot(n/norm(n), sum(A_v, 2))/2);

end