function [out] = q_pick_near(q, R)

    if q*R' < q*(-R)'
        out = -R;
    else
        out = R;
    end