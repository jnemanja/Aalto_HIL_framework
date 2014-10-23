
if exist('RunFile','var')==0
    clc
    close all
    clear
    WCDTorque
end

if exist('RunFile','var')==1

    % dimensions of satellite body
    x = [SatWidth 0 0]';              % dimension of satellite body x-axis in SBRF (x,y,z) [m]
    y = [0 SatLength 0]';              % dimension of satellite body y-axis in SBRF (x,y,z) [m]
    z = [0 0 SatHeigth]';              % dimension of satellite body z-axis in SBRF (x,y,z) [m]
    
    % arrays used to store max torques at each CoM
    m_torque_r = [];
    torque_r = [];

    for A=0:1:10
        for B=0:1:10
            for C=0:1:10

                % calculate the radiation direction unit vector
                v_hat = [A,B,C]';
                v_hat = v_hat/norm(v_hat);
                
                % calculate the projection matrix
                P = 1/(A^2 + B^2 + C^2)*[B^2+C^2, -A*B, -A*C;-A*B,A^2+C^2,-B*C;-A*C,-B*C A^2+B^2];
                
                % project the satellite body and solar array dimension onto the
                % plane ortogonal to the radiation direction unit vector
                xp = P*x;
                yp = P*y;
                zp = P*z;
                
                A_sat = norm(cross(xp,yp))+norm(cross(yp,zp))+norm(cross(xp,zp));
                F_r = Crk*A_sat*Prmf;
                
                torque_r = [torque_r,norm(cross(GoM-CoM,F_r*v_hat))];
            end
        end
    end
    
    % maximum torque is found
    m_torque_r = max(torque_r);
    
    disp('Max Radiation Torque')
    max(m_torque_r)
    
end
