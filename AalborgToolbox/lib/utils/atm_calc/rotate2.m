function rp = rotate2(pol, r)

    rp = angle2dcm(r(1)*3.14/180, r(2)*3.14/180, r(3)*3.14/180, 'XYZ')*pol;

end