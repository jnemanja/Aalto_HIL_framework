function y = PD(u)

    q = u(1:3);
    w = u(5:7);
    B = u(11:13);
    
    P = 5*10^-4;
    D = 2*10^-1;
    T = (-P*q-D*w);
    y = [[1 1 1]'.*cross(B,T)/(B'*B); 0; 0; 0; 0; 0; 0];

end

