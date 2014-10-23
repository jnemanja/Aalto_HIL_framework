function rp = rotate(pol, r)
    a = r(1)*3.14/180;
    b = r(2)*3.14/180;
    c = r(3)*3.14/180;
    R = [cos(b)*cos(c) -cos(b)*sin(c) sin(b);
         cos(a)*sin(c)+cos(c)*sin(a)*sin(b) cos(a)*cos(c)-sin(a)*sin(b)*sin(c)...
         -cos(b)*sin(a);
         sin(a)*sin(c)-cos(a)*cos(c)*sin(b) cos(c)*sin(a)+cos(a)*sin(b)*sin(c)...
         cos(a)*cos(b)];
    rp = R'*pol;
end