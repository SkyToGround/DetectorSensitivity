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
As quite a few arguments are required in order for the application to run, it is suggested that you start by testing the provided example-files. Run an example using the "--config_file" argument, e.g. "./MobileDetectorSim --config_file example_3.ini". Use the `--verbose` flag to get information on the progress of the calculations.

Note that the file "example_1.ini" has a set of settings which requires a significant amount of time to run.

##Disclaimer
I have run several tests of the code to make sure that the calculations are correct. However, I can not ensure that they really are and thus the results should be dubble checked with experimental data if it is to be used with real measurement systems.

Due to some recent changes in the code, when setting the application to find the optimal integration time, it is possible for the minima finding algorithm to get stuck on a local minima and thus the `act_vs_time` command should be used so that you can verify that the code actually found the global minima.

If "unexpected" input parameters are used, the code might very well enter an infinite loop. If you are running in `verbose` mode and have not received a new message in a while, the code has probably entered an infinite loop. The suggested solution to that (for now) is to test other input values or to start debugging.

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
* Check if it is possible to remove the use of `long double` in the list mode part of the code.
