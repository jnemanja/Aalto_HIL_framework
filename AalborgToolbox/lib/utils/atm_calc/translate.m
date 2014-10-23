function tp = translate(pol, v) 
    tp = pol+v'*ones(1, size(pol, 2));
end