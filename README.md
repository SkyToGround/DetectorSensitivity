#DetectorSim

This application is used to calculate the radioactive source activity which can be detected by a detector that is moving past it (or a source moving past a detector). A more through explanation of the function of this application can be found in the manual.

##Requirements
The following libraries are needed to compile the code:
* Boost
* GSL
* Eigen3

In addition, CMake (version >= 3.1.0) is required to create the make-file.

##Compiling
* Change directory to the source code directory (which contains the "CMakeLists.txt" file).
* Run CMake:
	CMake .
* Then compile the application:
	make
* When the application is compiled, you will find an executable in the current directory.
