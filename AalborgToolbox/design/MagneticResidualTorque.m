
if exist('RunFile','var')==0
    clc
    close all
    clear
    WCDTorque
end

%% Init
if exist('RunFile','var')==1

    % torque from a wire loop with dimensions 5.7 x 5.7 cm and a current of
    % 1[W]/3.3[V]. -> gives the same max magnetic residual torque 
    N_mag_res_imaginary=(1*1/3.3*0.057*0.057*maxB)

    Nmagsc_max=0.001*maxB; % 0.001 A*m^2 is 100 times less than what Wertz1 suggests.
    Nmagsc_min=0.001*minB;
    
    disp('Max Magnetic Residual (High Field Strength)')
    Nmagsc_max
    
    disp('Max Magnetic Residual (Low Field Strength)')
    Nmagsc_min
end

