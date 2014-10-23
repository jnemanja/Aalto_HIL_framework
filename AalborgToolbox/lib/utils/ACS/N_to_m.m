function [ output ] = N_to_m(input)
% This function Projects the ideal control torque down to the power optimal
% magnetic moment, which lies in the plane perpendicular to the geomagnetic
% field
%
% Written by Kasper V. and Kasper F. Jensen (AAU - 2010)

persistent msg;

if isempty(msg)
    msg=[];
end

N_S=input(1:3);
B_S=input(4:6);

if max(isnan(B_S)) == 1
    disp('B_S was NaN');
    msg=[msg,'B_s was NaN : '];
    assignin('base','msg',msg);
    m_ctrl=[0;0;0];
else
    B_S_skew=[0        -B_S(3)  B_S(2);
              B_S(3)   0        -B_S(1);
              -B_S(2)  B_S(1)   0     ];
    
    % This is the projection proposed by Rafael Wisniewski 
    m_ctrl=(B_S_skew*N_S/(B_S(1)^2+B_S(2)^2+B_S(3)^2));
end

%prik1=(N_S./norm(N_S))'*(B_S./norm(B_S));
%angle_NS_to_BS=acos(prik1)

%N_S_mag=B_S_skew'*m_ctrl;
%prik2=(N_S_mag./norm(N_S_mag))'*(B_S./norm(B_S));
%angle_NSmag_to_BS=acos(prik2)

%norm1=norm(N_S)
%norm2=norm(N_S_mag)
%ratio=norm1/norm2

output=m_ctrl;

end

