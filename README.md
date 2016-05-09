#DetectorSim
This application is used to calculate the radioactive source activity which can be detected by a detector that is moving past it (or a source moving past a detector). A more through explanation of the function of this application can be found in the manual (which is not done yet). Note that I leave you no guarantee that the application will work as expected. Use it at your own risk.

##Requirements
The following libraries are needed to compile the code:
* Boost
* GSL
* Eigen3

In addition, CMake (version >= 3.1.0) is required to create the make-file.

##Compiling
Start by chaning directory to the source code directory (which contains the "CMakeLists.txt" file).
* Run CMake: `CMake .`
* Then compile the application: `make`

When the application is compiled, you will find an executable in the current directory.

##Running the examples
As quite a few arguments are required in order for the application to run, it is suggested that you start by testing the provided example-files. Run an example using the "--config_file" argument, e.g. "./MobileDetectorSim --config_file example_3.ini".

Note that the file "example_1.ini" has a set of settings which requires a significant amount of time to run.


##To-Do
These are in no particular order.
* Fix problem related to the algorithm finding a local minima instead of a global one.
* Add simulation of detectors from -90 to 90 degrees.
* Improve verbose mode.
* Optimise the minima finding algorithm when using list mode (maybe).
* Write a manual.
* Take 3 axes into account when calculating distance (to make useful for calculations on airplanes).
* Take the attenuation of gamma radiation in air into account.
* Pick the initial C_L value with more intelligence. Too low values of C_L increases the calculation time quite a bit.
* Better integration algorithm.
