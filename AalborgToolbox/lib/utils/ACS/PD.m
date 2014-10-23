function y = PD(u)

    q = u(1:3);
    w = u(5:7);

    P = 0*10^-4;
    D = 2*10^-4;
    
    y = - P*q - D*w;

end

