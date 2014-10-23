function startup(MEX)
%this script is automatically run when matlab is started in this dir
% it adds the following to the matlab path
sub_dirs = {
    'lib/'
    'lib/AlbedoToolbox-1.0'
    'lib/AlbedoToolbox-1.0/refl_data'
    'lib/AlbedoToolbox-1.0/refl_data/2005'
    'lib/utils'
    'lib/utils/igrf'
    'lib/utils/atm_calc'
    'lib/utils/ephemeris'
    'lib/utils/math_utils'
    'lib/utils/ACS'
    'lib/pictures'
    'design/'
    'models/'
    'test/'
    'plots/'
    'mex_files/'
    ''
	   };

  source_files = {
    '../lib/utils/scdynamics.c'
    '../lib/utils/zonal.c'
    '../lib/utils/ephemeris/sgp4S.c'
    '../lib/utils/igrf/igrfS2.c'
    '../lib/utils/math_utils/qmul.c'
    '../lib/utils/math_utils/qestimator.c'
	   };
  
len=length(sub_dirs);
disp('Adding to path...')
for i= 1:len,
  temp = [pwd '/' char(sub_dirs(i))];
% disp(temp)
  addpath(temp)
end
disp('Done!')

clear sub_dirs len i temp
if nargin~=0
disp('Compiling mex files...')
% mkdir mex_files
cd mex_files
len=length(source_files);
for i= 1:len,
  temp = [pwd '/' char(source_files(i))];
  disp(['Compiling ' temp]);
  eval(['mex' ' ' temp]);
end
cd ..
disp('Done!')
end
