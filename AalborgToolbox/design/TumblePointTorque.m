
if exist('RunFile','var')==0
    clc
    close all
    clear
    WCDTorque
end

%% Init
if exist('RunFile','var')==1

    % Detumble torque  
    Ndet=eig(diag(Isat_S))*((InitDetumble-EndDetumble)/(norbits*Torbit));
    Ndet=max(max(Ndet))

    % Pointing torque
    phi=asin(sqrt((ERadiusMean+AltitudeSat)^2-ERadiusMean^2)/(ERadiusMean+AltitudeSat));
    Tpoint=Torbit/(2*pi)*phi %time for half a ground station pass

    Npoint=eig(diag(Isat_S))*(omegaPoint/Tpoint);

    Npoint=max(max(Npoint))*2 % Factor 2 to handle better than the requirement
end