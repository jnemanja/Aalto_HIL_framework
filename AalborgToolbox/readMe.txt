%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% AAUSAT3 Simulator  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Installation Guide:

Requires:
    - SeDuMi (http://sedumi.ie.lehigh.edu/).
    - YALMIP (http://users.isy.liu.se/johanl/yalmip/).
    - AAUSAT3 library files.
    - Matlab mex compatible compiler:
	- Visual Studio 2008 Professional Edition [TESTED] (http://www.microsoft.com/downloads/details.aspx?displaylang=en&FamilyID=83c3a1ec-ed72-4a79-8961-25635db0192b)
	- GCC compiler (linux).
	- Microsoft Windows Software Development Kit (http://www.microsoft.com/downloads/) [Windows only].

Steps (windows):
1. Install Visual Studio 2008 Professional Edition.
2. Install Microsoft Windows SDK.
3. Setup the mex compiler in Matlab using ''mex -setup'' command.
4. Setup the Matlab compiler in Matlab using ''mbuild -setup'' command.
5. Install YALMIP.
6. Install SeDuMi.
7. Remember to add paths to Matlab.
8. Run the ''install_sedumi'' command from the SeDuMI folder.
9. Test whether YALMIP is working using 'yalmiptest'.

Steps (linux):
1. Setup the mex compiler in Matlab using ''mex -setup'' command.
2. Setup the Matlab compiler in Matlab using ''mbuild -setup'' command.
3. Install YALMIP.
4. Install SeDuMi.
5. Remember to add paths to Matlab.
6. Run the ''install_sedumi'' command from the SeDuMI folder.
7. Test whether YALMIP is working using 'yalmiptest'.

Initialize environment:
1. Go to the AAUSAT3 folder.
2. Run the ''startup command'' (to recompile mex-files run ''startup('mex')'').

-For automatic setup of path, start matlab in the path of startup.m, which sets up the paths.

The simulator source files are sorted in designated directories, which name reflects the purpose
for which they are used. The different directories are 
#-------------------------#-----------------------------------------------------------------------------#
|	Directory Name 	  |	Content Description							              |
#-------------------------#-----------------------------------------------------------------------------#
|	root  		  |	Files used to ensure simulation on various versions of Matlab, 		  |
|			        |	the magnetic field model data and startup.m      			        |
#-------------------------#-----------------------------------------------------------------------------#
|	lib		        |	Contains a simulink library with the components of the simulator.	        |
#-------------------------#-----------------------------------------------------------------------------#
|	lib/albedo_toolbox  | 	Contains a source model for the albedo toolbox used in			  |
|      			  |	the modeling developed by M.Sc.E.E., Ph.D., Dan Bhanderi, AAU.		  |
#-------------------------#-----------------------------------------------------------------------------#
|	lib/pictures  	  |	Contains jpeg-files used in different simulink block masks. 		  |
#-------------------------#-----------------------------------------------------------------------------#
|	lib/utils   	  |	Contains mex- and m-files used in the mathematical modeling of		  |
|     			  |	the satellite sorted into subdirectories.				              |
#-------------------------------------------------------------------------------------------------------#
|	misc 		        |	Miscellaneous script files for, e.g , compiling mex-files		        |
|     			  |	printing simulink models or illustrating vectors.			        |
#-------------------------#-----------------------------------------------------------------------------#
|	design 		  |	Files used to create different actuators, estimators and controllers.	  |
#-------------------------#-----------------------------------------------------------------------------#
|	models 		  |	Contains various simulink model files used as templates for test cases.   |
#-------------------------#-----------------------------------------------------------------------------#
|	test 		        |	Contains simulink model files with specific controllers			  |
|			        |	used to evaluate the individual controllers.				        |
#-------------------------#-----------------------------------------------------------------------------#

Generally a description of the m-files used in the simulator can also be found in the individual files.