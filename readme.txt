 This folder contains the source file as well as the build project under different systems with different settings

Requirements:
To compile PAWS, you need 
(a) Matlab installed. -- this is a condition we are trying to get rid of in the future.

To run PAWS, you need
either Matlab, or Matlab Compiler Runtime (freely available), see below for set up
(the netcdf dependencies have been succesfully removed. Goodbye)
 

/src: source files
   ../PAWS: PAWS alone source files
   ../CLM: original and modified CLM source files
   ../include: include directory containing header files for compilation. preproc.h can be edited to suit the purpose of the build. But other files should not be modified.
   ../include_PAWS_only: If only compiling PAWS with simple vegetation module, without CLM functionalities, then this is the folder to be put under 'Additional Include Directories' (or -I command option)
   ../bin: two files that come with nedcdf, not even sure what they are for.
   ../obsolete: outdated source files which may occasionally be used during testing stage
   ../lib: some library files that are needed during compilation. On Windows, netcdf.dll need to be copied to a system paths (e.g. C:\Windows\system32). On linux, the netcdf paths needs to be known to bash at runtime.
/PAWS: Microsoft Visual Studio project files of PAWS together with CLM. Inside the folder there are debug and release build. The project file has all paths,include/library directories already defined.
/PAWS_FORTRAN: a Microsoft Visual Studio build with only PAWS functionalities, with PAWS's own simple vegetation module. no CLM functionalities.

If your matlab root is not:
C:\Program Files\MATLAB\R2009a
then to use the existing build project files, your $MATLABROOT must be correctly updated in 
Project Properties --> Fortran --> Preprocessing --> Additional Include Directories
Project Properties --> Linker  --> General       --> Additional Library Directories

Test git
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I've been using Win64 version of the build under PAWS_debug. The Win32 version may be slightly outdated. A couple of files may need to be added into "Source Files". Other settings should be working
contact me if you need support on Win32 cshen@lbl.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BUILD PAWS on Windows:
/include directory needs to be added into Additional Include Directories. If multiple directories are needed, use semi-colon to separate them, there must be no blanks (same applies to Additional Library Directories) after colons.
Also add:
C:\Program Files\MATLAB\R2009a\extern\include;

(previously used library files have now been removed as we removed netcdf dependencie)
Also add:
C:\Program Files\MATLAB\R2009a\extern\lib\win32\microsoft
(on Win64: C:\Program Files\MATLAB\R2010b\extern\lib\win64\microsoft)
Add these into Input/Additional Dependencies even if they are not in the /lib folder:
libmat.lib libmx.lib

Other Compiler Options:
Preprocess Source File: Yes (/fpp)        [Fortran -> Preprocessor]
Check Routine Interfaces: No              [Fortran -> Diagnostics] (there should be no (/warn:interfaces) on the command line)
Disable Default Library Search Rules: No  [Fortran -> libraries]
Check Array and String Bounds:        No  [Fortran -> Run-time]
Runtime Library: Debug Multithreaded or Mutlithreaded

LINKER, Command Line: /stack:20000000





=========================
== BELOW IS NOT CURRENTLY RELEVANT
we can compile all programs into a static/dynamic library, and at the last step mex it with the gateway subroutine

ifort VData.F90 SnowLSub.f snow.F90 Flow.F90 /c /MD 

(the /MD flag is to use the multi-threaded library as Matlab. see error 4098
http://www.mathworks.com/support/solutions/en/data/1-12L4PQ/?solution=1-12L4PQ
http://msdn.microsoft.com/en-us/library/Aa267384
)

Then build a library with the object files:
LIB /out:flow.lib VData.obj snowLSub.obj snow.obj flow.obj

this will produce a library file called flow.lib, last step is to link it with mex file:

mex prism.F90 flow.lib -output prism

In this manner, we are no longer confined by the fortran source file F77 conventions. We can compile all relevant codes into a library that is linked.
