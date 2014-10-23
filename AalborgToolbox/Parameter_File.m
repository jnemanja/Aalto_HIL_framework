% Constants and Parameters for ESTCube-1 Simulation
% =================================================

% This file contains all parameters needed to simulate the detumbling phase
% of ESTCube-1. The Simulation bases on the Simulink(R)-Testenviroment,
% developed by the AAUSAT-Team from Aalborg, Denmark.
% Because of the differences between the satellites some changes on the
% model have been done.


% Index of Contents
% =================

% 1 ESTCube-1 Satellite Parameters
% 2 Magnettorquers & Drivers Parameters
% 3 Measurements Parameters
% 4 TLE
% 5 B-Dot Contoller Parameters
% 6 Controller
% 7 Changelog


% 1 ESTCube-1 Satellite Parameters
% ================================

Center_of_Mass_sat=[0.05,0.05,0.055];                            % [m]
% Center of Mass based on estimations (geometrical center)
% 06.04.2011 - Final data from STR (SolidWorks-Model)

initial_Temperature_sat=293.15;                                  % [K]
% Temperature when detumbling phase starts - source AAUSAT-III
% 06.04.2011 - TBD

Mass_sat=1.33;                                                   % [kg]
% Max. acceptable mass (CubeSat-Standarts)
% 06.04.2011 - Final data needs to be collected from all Subsystems

Dimension_sat=[0.1,0.1,0.11];                                    % [m]
% Proposed dimensions - CubeSat-Standarts
% 06.04.2011 - TBC

Inertia_Matrix_sat=[0.002,0.0022,0.004];                         % [kg*m^2]
% Estimations basing on early calculations (Phase A)
% 06.04.2011 - Final data from STR (SolidWorks-Model)

Initial_Time=2455153.1954318;                                   % [JD Days]
% Add 0.0169 to satellite over north pole
% Initial time as Julian Date
% Current data from AAUSAT-III
% 07.04.2011 - TBC

Controller_Reference_Frame=[0 0 0 1];                            % [1]
% Controller Reference Frame in quaternions
% Controller Reference Frame = Spacecraft Reference Frame
% 27.09.2011 - TBC
%0 0 0.1736 0.9848
Initial_Attitude = [0 0 0 1];                                      % [1]
% Initial Attitude in quaternions
% 07.04.2011 - TBD/TBC
% [0 0 30] = [0 0 -0.2588 0.9659]
% [30 0 0] = [-0.2588 0 0 0.9659]
% [0 30 0] = [0 -0.2588 0 0.9659]
% [89 0 0] = [-0.7009 0 0 0.7133]
% [0 89 0] = [0 0.7009 0 0.7133]
% [90 0 0] = [-0.7071 0 0 0.7071]
% [0 0 -90] = [0 0 -0.07071 -0.7071]
% [180 0 0] = [-1 0 0 0]

Initial_Angular_Rate=[0,0,0]*pi/180;                          % [rad/s]
% Values can be changed for different simulation scenarios
% 07.04.2011 - TBC

Atmospheric_Drag_Coefficient=2;                                  % [1]
% Source AAUSAT-III
% 07.04.2011 - TBC

Atmospheric_Density=5e-14;                                       % [kg/m^3]
% Source AAUSAT-III
% 07.04.2011 - TBC

Solar_Momentum_Flux=4.5565e-6;                                   % [kg/ms^2]
% Source AAUSAT-III
% 07.04.2011 - TBC

Solar_Absorbtion_Coefficient=1.5;                                % [1]
% Source AAUSAT-III
% 07.04.2011 - TBC

Max_Angle_of_Theta=65;                                           % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBD/TBC

Reduction_Factor=1;                                              % [1]
% Source AAUSAT-III
% 07.04.2011 - TBC

Albedo_Samplin_Time=1;                                           % [s]
% Source AAUSAT-III
% 07.04.2011 - TBC


% 2 Magnettorquers & Drivers Parameters
% =====================================

voltage_coil_bang_bang=5;                                        % [V]
% Coil Voltage - just used for the "Bang-Bang-Controller"
% 24.09.2011

Coil_voltage = 4;
% Battery typical voltage. 11.12.2011                            % [V]

Dimensions_coil=[0.078,0.057];                                   % [m]
% Coil average dimensions [x,y]. EC1-C-ADCS-V1-Coil_Calculator. 11.12.2011
% Dimensions_coil=[0.08,0.08]; 07.04.2011. 

Number_of_Windings_coil=400;                                     % [1]
% EC1-C-ADCS-V1-Coil_Calculator. 11.12.2011

Wire_Cross_Sectional_Area_coil=0.0000000283;                     % [m^2]
% EC1-C-ADCS-V1-Coil_Calculator. Wire diameter 19mm. 11.12.2011.

Magnetic_Momentum_coil=0.116;                                    % [A*m^2]
% EC1-C-ADCS-V1-Coil_Calculator. 4V, 57Ohm, [78,57] mm. 11.12.2011.

Wire_Resistivity_coil=1.68e-8;                                   % [Ohm*m]
% Specific Resistivity of copper - based on
% http://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity
% 07.04.2011 - TBC

Wire_Resistivity_Temperature_Coefficient_coil=3.9e-3;            % [(deg C)^-1]
% Wire Resistivity Temperature Coefficient of copper - Source AAUSAT-III
% 07.04.2011 - TBC

Wire_Resistivity_Base_Temperature_coil=293.15;                   % [K]
% Source AAUSAT-III
% 07.04.2011 - TBC

Actuator_Frequency_coil=10;                                      % [Hz]
% Needs to be alligned to "Magnetometer_Sampling_Frequency_sens"
% Source AAUSAT-III
% 07.04.2011 - TBD

Actuator_Duty_Cycle_coil=88;                                     % [%]
% Source AAUSAT-III
% 07.04.2011 - TBD


% 3 Measurements Parameters
% =========================

Magnetometer_Sampling_Frequency_sens=10;                         % [Hz]
% Needs to be alligned to "Actuator_Frequency_coil"
% Source AAUSAT-III
% 07.04.2011 - TBD

Magnetometer_Noise_sens=3;                                       % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBC

Magnetometer_Bias_sens=[0 0 0];                                  % [nT]
% Source AAUSAT-III
% 07.04.2011 - TBC

Magnetometer_Placement_Error_heading_sens=0;                     % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBC

Magnetometer_Placement_Error_tilt_sens=0;                        % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBC

Sun_Sensor_Sampling_Frequency_sens=1;                            % [Hz]
% Source AAUSAT-III
% 07.04.2011 - TBD

Sun_Sensor_Noise_sens=3.33;                                      % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBC

Sun_Sensor_Bias_sens=[0 0];                                      % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBC

Number_fo_Sun_Samples_sens=1;                                    % [1]
% Source AAUSAT-III
% 07.04.2011 - TBD

Gyroscope_Sampling_Frequency_sens=1;                             % [Hz]
% Source AAUSAT-III
% 07.04.2011 - TBD

Gyroscope_Noise_sens=1;                                          % [deg/s]
% Source AAUSAT-III
% 07.04.2011 - TBC

Gyroscope_Bias_sens=[0 0 0];                                     % [deg/s]
% Source AAUSAT-III
% 07.04.2011 - TBC

Gyroscope_Placement_Error_heading_sens=0;                        % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBC

Gyroscope_Placement_Error_tilt_sens=0;                           % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBC


% 4 TLE
% =====

Epoch_tle=09321.69542252;                                        % [YYDDD.DDDDDDDD]
% Source AAUSAT-III
% 07.04.2011 - TBD

Mean_Motion_1st_tle=0.00000555;                                  % [rev/day^2]
% Source AAUSAT-III
% 07.04.2011 - TBD

Mean_Motion_2nd_tle=0.0;                                         % [rev/day^3]
% Source AAUSAT-III
% 07.04.2011 - TBD

B_Drag_tle=.76898e-4;                                            % [Eart radii^-1]
% Source AAUSAT-III
% 07.04.2011 - TBD

Inclination_tle=97.9251;                                         % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBD 

Right_Ascention_of_Ascending_Node_tle=24.1597;                   % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBD

Eccentricity_tle=0.0014245;                                      % [1]
% Source AAUSAT-III
% 07.04.2011 - TBD

Argument_of_Perigee_tle=254.1774;                                % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBD

Mean_Anomaly_tle=105.7862;                                       % [deg]
% Source AAUSAT-III
% 07.04.2011 - TBD

Mean_Motion_tle=14.81823818;                                     % [rev/day]
% Source AAUSAT-III
% 07.04.2011 - TBD

Sampling_Time_tle=1;                                             % [sec]
% Source AAUSAT-III
% 07.04.2011 - TBC

% 4.A Other orbital parameters
% ============================
Earth_Mass = 5.9722e+024;                                        % [kg]
Gravitational_Constant = 6.6738e-011;                            % [m^3*kg^(-1)*s(-2)]
Earth_Radius = 6371000;                                          % [m]
Mu = Gravitational_Constant * Earth_Mass;                        % [m^3*s(-2)] Standard gravitational parameter
Day_Length = 86400;                                              % [s]
Semimajor_Axis=((Day_Length/(Mean_Motion_tle*2*pi))^2*Mu)^(1/3); % [m]
% General formula is Semimajor_Axis=(Mu/(n^2))^(1/3) where n is Mean Motion
% but since in therory Mean Motion is mesured in rad/sec then the
% transformation from rev/dey is needed.
Orbital_Period = 2*pi*sqrt(Semimajor_Axis^3/(Mu));               % [s]
% Formulas are taken from Fundamentals of Astrodynamics and Applications,
% David A. Vallado

% 5 B-Dot Contoller Parameters
% ============================

Gain_b_dot=25000;                                                 % [1]
% Controller gain
% 23.09.2011 - TBD

Filter_cutoff_frequency_b_dot=0.7;                               % [1]
% Break frequence - source AAUSAT-III
% 07.04.2011 - TBD

Sample_time_b_dot=0.1;                                           % [s]
% Source AAUSAT-III
% 07.04.2011 - TBD


% 6 Controller
% ============

K_initial=[-18.5,0,0,-14500,0,7750;
            0,-32,0,0,-9500,0;
            0,0,-14.5,2200,0,-16000];
% Try and Error
% 27.09.2011 - TBD

K_nadir=[-18.5,0,0,-14500,0,7750;
            0,-32,0,0,-9500,0;
            0,0,-14.5,2200,0,-16000];
% Try and Error
% 27.09.2011 - TBD

K_spinup=[-18.5,0,0,-14500,0,7750;
            0,-32,0,0,-9500,0;
            0,0,-14.5,2200,0,-16000];
% Try and Error
% 27.09.2011 - TBD

% 7 Spin-Up Controller from the Article
% ============

W=[1,0,0; 0,1,0; 0,0,1];

k=0.004;
%k=0.04
k_1=1.5;
%k_1 = 2;
k_2=0.5;
%k_2 = 0.25

% P for JC2Sat
% P=[1,0,0; 0,0,0; 0,0,1];

% P for ESTCube-1
P=[1,0,0; 1,0,0; 0,0,0];

% Spin axis
spin_axis = [0; 0; 1];
spin_axis_backup = [0; 1; 0];

%ESTCube-1
Inertia_Matrix_sat2=[0.002, 0, 0; 0, 0.0022, 0; 0, 0, 0.004];  

%JC2Sat
%Inertia_Matrix_sat2=[1.0669, -0.0258, 0.0000; -0.0258, 1.0322, 0.0033; 0.0000, 0.0033, 1.0475];                         % [kg*m^2]

%JC2Sat inertia matrix in principal coordinates
%Inertia_Matrix_sat2=[1.0178, 0, 0; 0, 1.0792, 0; 0, 0, 1.0477];                         % [kg*m^2]


% Spin-Up article
% 5.12.2011

omega_d_scal = 360*pi/180;     % [rad/s]
omega_d = spin_axis*omega_d_scal;
%omega_d = [0;0;omega_d_scal]; % [rad/s]
%omega_d = [360;0;0]*pi/180;  % [rad/s]
% 1 rotation per second

m_mt_max = 0.127;              % [A*m^2]
% Old value for magnetotorquers


% 8 Changelog
% ===========

% 06.04.2011 - Tobias Scheffler
% - designations updated
% - some descriptions added

% 07.04.2011 - Tobias Scheffler
% - descriptions added
% - Voltage_coil added
% - TLE added

% 23.09.2011 - Tobias Scheffler
% - Voltage_coil removed
% - Gain_b_dot added

% 24.09.2011 - Tobias Scheffler
% - voltage_coil_bang_bang added

% 5.12.2011 - Erik Kulu
% Spin-Up parameters from article and elsewhere

% 11.12.2011 - Erik Kuku
% Magnetorquer parameters changed