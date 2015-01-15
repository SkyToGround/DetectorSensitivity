//
//  main.cpp

#include <iostream>
#include <eigen3/Eigen/Dense>
#include "Detector.h"
#include <string>
#include <vector>
#include <eigen3/unsupported/Eigen/FFT>
#include <boost/timer/timer.hpp>
#include <Python.h>
#include "Response.h"
#include "Measurement.h"
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/shared_ptr.hpp>
#include "CalcCoordinator.h"
#include <boost/program_options.hpp>

#include <boost/timer/timer.hpp> //Fix me; remove

namespace po = boost::program_options;

using namespace std;

void plotData(string plotFunc, vector<double> xvec, vector<double> yvec);
void plotData(string plotFunc, const ArrayXd &x, const ArrayXd &y);
void CalcActTest();
void MathExplanationPlots();
void TimeVsActExamplePlot();
void PerformCalculations();
//int main2(int argc, const char * argv[]);

void limitedCalcTest();


int main2(int argc, const char * argv[])
{
	//TimeVsActExamplePlot();
	MathExplanationPlots();
	//PerformCalculations();
	//limitedCalcTest();
	return 0;
}

int main(int argc, const char * argv[])
{
	po::options_description desc("Possible arguments");
	desc.add_options()
	("help", "Show this help message")
	("fph", po::value<double>()->default_value(1.0), "Acceptable number of false positives per hour. Must be > 0. Default value is 1.")
	("beta", po::value<double>()->default_value(0.05), "Probability of false negative. Must be a value > 0 and < 1. Default value is 0.05.")
	("uncertainty", po::value<bool>()->default_value(false), "Perform uncertainty calculation. Possible values are 'true' and 'false'. Default is 'false'.")
	("uncertainty_loops", po::value<unsigned int>()->default_value(100), "Number of loops in uncertainty calculation. Default value is 100.")
	("calc_type", po::value<string>(), "Type of calculation to perform. Possible values are: 'mean', 'best' and 'worst'. Accepts multiple options. Defaults to 'mean'.")
	("output", po::value<string>(), "If argument is given, will output to given file in JSON-format. If not; outputs results to screen.")
	("distance", po::value<double>(), "The distance at which the source is placed in metres. Must be > 0. Accepts multiple values. Default value is 10.")
	("velocity", po::value<double>(), "The velocity of the detector on metres per second. Must be > 0. Accepts multiple values. Default value is 8.33 (30 kph).")
	
	("dist_response", po::value<string>(), "Detector response input. Must in the following format: \"distance:pulses:time\" Where \"distance\" is the distance between source and detector in metres, \"pulses\" is the number of pulses reigstered by the detector and \"time\" is the live time of the measurement. Accepts multiple values. See manual and examples for more information.")
	
	("ang_response", po::value<string>(), "Detector response input. Must in the following format: \"angle:pulses:time\" Where \"angle\" is the ????, \"pulses\" is the number of pulses reigstered by the detector and \"time\" is the live time of the measurement. Accepts multiple values. See manual and examples for more information.")
	;
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if (vm.count("help") or argc == 1) {
		cout << desc << "\n";
		return 1;
	}
	
	return 0;
}

void TimeVsActExamplePlot() {
//	vector<string> dists{"5", "10", "15"};
//	Detector det10Am("pvc10cm_bkg", "pvc10cm_avst_@m_am", dists, "pvc10cm_@_10m_am", 1.1, 0.11, false);
//	
//	det10Am.SetDistance(10.0);
//	int nrOfTimes = 2000;
//	ArrayXd testTimes = ArrayXd::LinSpaced(nrOfTimes, 1.0, 15);
//	ArrayXd resM(nrOfTimes);
//	vector<Detector::CalcType> calcTypes{Detector::CalcType::BEST, Detector::CalcType::MEAN, Detector::CalcType::WORST};
//	for (int w = 0; w < testTimes.size(); w++) {
//		double alpha = testTimes[w] / 3600;
//		det10Am.SetIntegrationTime(testTimes[w]);
//		double c_M = det10Am.CalcActivity(alpha, 0.05, calcTypes[1]);
//		resM[w] = c_M;
//	}
//	int M_loc;
//	double bestM = resM.minCoeff(&M_loc);
//	double bestTime = testTimes[M_loc];
//	cout << "Best activity is " << bestM << " which was found at time " << bestTime << " s." << endl;
//	plotData("ExamplePlot", testTimes, resM);
	
	vector<string> dists{"5", "10", "15"};
	//Detector det("he3_bkg", "he3_am_@m_avst", dists, "he3_am_10m_@", 1.1, 0.11, true);
	
	Detector det("pvc10cm_bkg", "pvc10cm_avst_@m_am", dists, "pvc10cm_@_10m_am", 1.1, 0.11, false);
	
	det.SetDistance(10.0);
	int nrOfTimes = 2000;
	ArrayXd testTimes = ArrayXd::LinSpaced(nrOfTimes, 1.0, 15);
	ArrayXd resM(nrOfTimes);
	vector<Detector::CalcType> calcTypes{Detector::CalcType::BEST, Detector::CalcType::MEAN, Detector::CalcType::WORST};
	for (int w = 0; w < testTimes.size(); w++) {
		double alpha = testTimes[w] / 3600;
		det.SetIntegrationTime(testTimes[w]);
		double c_M = det.CalcActivity(alpha, 0.05, calcTypes[1]);
		resM[w] = c_M;
	}
	int M_loc;
	double bestM = resM.minCoeff(&M_loc);
	double bestTime = testTimes[M_loc];
	cout << "Best activity is " << bestM << " which was found at time " << bestTime << " s." << endl;
	plotData("ExamplePlot", testTimes, resM);
}

void MathExplanationPlots() {
	vector<string> distsPu{"2.5", "5", "10"};
	Detector detPu("he3_bkg", "he3_pu_@m_avst", distsPu, "he3_pu_5m_@", 0.1, 0.05, true);
	
	detPu.SetIntegrationTime(1.0);
	detPu.SetDistance(2.0);
	detPu.SetVelocity(4.0);
	
	ArrayXd times = ArrayXd::LinSpaced(1000, -3.0, 3.0);
	ArrayXd signal = detPu.S(times);
	//plotData("PythonPlotOld3", times, signal);
	plotData("TimeDescPlot", times, signal);
	return;
	detPu.SetDistance(2.0);
	int nrOfTimes = 2000;
	ArrayXd testTimes = ArrayXd::LinSpaced(nrOfTimes, 1.0, 7);
	ArrayXd resM(nrOfTimes);
	vector<Detector::CalcType> calcTypes{Detector::CalcType::BEST, Detector::CalcType::MEAN, Detector::CalcType::WORST};
	for (int w = 0; w < testTimes.size(); w++) {
		double alpha = testTimes[w] / 3600;
		detPu.SetIntegrationTime(testTimes[w]);
		double c_M = detPu.CalcActivity(alpha, 0.05, calcTypes[1]);
		resM[w] = c_M;
	}
	int M_loc;
	double bestM = resM.minCoeff(&M_loc);
	double bestTime = testTimes[M_loc];
	cout << "Best activity is " << bestM << " which was found at time " << bestTime << " s." << endl;
	double timeTemp = 1.0, mTemp;
	for (int i = 0; i < 10; i++) {
		FindLocalMinima(detPu, 1.0, 0.05, calcTypes[1], timeTemp, mTemp);
		cout << "Found minima at t = " << timeTemp << ", where M = " << mTemp << endl;
		timeTemp = timeTemp + 0.1;
	}
	
	plotData("MinimisePlot", testTimes, resM);
}

void limitedCalcTest() {
	CalcCoordinator calc(OutputType::SCREEN, string("results.txt"));
	vector<string> dists{"5", "10", "15"};
	
	Detector detHeAm("he3_bkg", "he3_am_@m_avst", dists, "he3_am_10m_@", 1.1, 0.11, true);
	
	detHeAm.SetDistance(10.0);
	detHeAm.SetIntegrationTime(1.0);
	detHeAm.SetVelocity(8.333);
	
	cout << "Ref. signal: " << detHeAm.S(2.0) << endl;
	
	vector<double> calcDists{2.5, 5.0, 10.0, 20.0};
	vector<string> calcNames{"best", "mean", "worst"};
	
	string dist_string = string("d") + lexical_cast<string>(int(calcDists[1] * 100)) + string("cm.");
	for (int i = 0; i < 1; i++) {
		calc.AddCalculation(detHeAm, 5.0, Detector::CalcType::MEAN, string("he3.am.") + dist_string + calcNames[1] + string("."), 1.0, 0.05, 0);
	}
	{
		boost::timer::auto_cpu_timer t;
		calc.RunCalculations();
		//10.343482s wall, 8.480000s user + 0.090000s system = 8.570000s CPU (82.9%) <- original
	}
	
}

void PerformCalculations() {
	CalcCoordinator calc(OutputType::JSON_FILE, string("results.txt"));
	vector<string> distsPu{"2.5", "5", "10"};
	vector<string> dists{"5", "10", "15"};
	
	Detector detHePu("he3_bkg", "he3_pu_@m_avst", distsPu, "he3_pu_5m_@", 4.74, 0.474, true);
	Detector detHeAm("he3_bkg", "he3_am_@m_avst", dists, "he3_am_10m_@", 1.1, 0.11, true);
	Detector detHeCf("he3_bkg", "he3_cf_@m_avst", dists, "he3_cf_10m_@", 6.63, 0.663, true);
	
	Detector det5Pu("pvc5cm_bkg", "pvc5cm_avst_@m_pu", distsPu, "pvc5cm_@_5m_pu", 0.1, 0.05, false);
	Detector det5Am("pvc5cm_bkg", "pvc5cm_avst_@m_am", dists, "pvc5cm_@_10m_am", 1.1, 0.11, false);
	Detector det5Cf("pvc5cm_bkg", "pvc5cm_avst_@m_cf", dists, "pvc5cm_@_10m_cf", 1.2, 0.36, false);
	
	Detector detBPu("bara_nai_bkg", "bara_nai_avst_@m_pu", distsPu, "bara_nai_@_5m_pu", 0.1, 0.05, false);
	Detector detBAm("bara_nai_bkg", "bara_nai_avst_@m_am", dists, "bara_nai_@_10m_am", 1.1, 0.11, false);
	Detector detBCf("bara_nai_bkg", "bara_nai_avst_@m_cf", dists, "bara_nai_@_10m_cf", 1.2, 0.36, false);
	
	Detector detUNPu("uppner_bkg", "uppner_avst_@m_pu", distsPu, "uppner_@_5m_pu", 0.1, 0.05, false);
	Detector detUNAm("uppner_bkg", "uppner_avst_@m_am", dists,"uppner_@_10m_am", 1.1, 0.11, false);
	Detector detUNCf("uppner_bkg", "uppner_avst_@m_cf", dists, "uppner_@_10m_cf", 1.2, 0.36, false);
	
	Detector det10Pu("pvc10cm_bkg", "pvc10cm_avst_@m_pu", distsPu, "pvc10cm_@_5m_pu", 0.1, 0.05, false);
	Detector det10Am("pvc10cm_bkg", "pvc10cm_avst_@m_am", dists, "pvc10cm_@_10m_am", 1.1, 0.11, false);
	Detector det10Cf("pvc10cm_bkg", "pvc10cm_avst_@m_cf", dists, "pvc10cm_@_10m_cf", 1.2, 0.36, false);
	
	
	
	vector<double> calcDists{2.5, 5.0, 10.0, 20.0};
	vector<string> calcNames{"mean", "best", "worst"};
	vector<Detector::CalcType> calcTypes{Detector::CalcType::MEAN, Detector::CalcType::BEST, Detector::CalcType::WORST};
	
	int uncertLoops = 100;
	double fpPerHour = 1.0;
	double beta = 0.05;
	
	string dist_string;
	for (int i = 0; i < calcDists.size(); i++) {
		for (int j = 0; j < calcTypes.size(); j++) {
			dist_string = string("d") + lexical_cast<string>(int(calcDists[i] * 100)) + string("cm.");
			calc.AddCalculation(detHePu, calcDists[i], calcTypes[j], string("he3.pu.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			calc.AddCalculation(detHeAm, calcDists[i], calcTypes[j], string("he3.am.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			calc.AddCalculation(detHeCf, calcDists[i], calcTypes[j], string("he3.cf.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			
			calc.AddCalculation(det5Pu, calcDists[i], calcTypes[j], string("pvc5cm.pu.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			calc.AddCalculation(det5Am, calcDists[i], calcTypes[j], string("pvc5cm.am.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			calc.AddCalculation(det5Cf, calcDists[i], calcTypes[j], string("pvc5cm.cf.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			
			calc.AddCalculation(detBPu, calcDists[i], calcTypes[j], string("onlynai.pu.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			calc.AddCalculation(detBAm, calcDists[i], calcTypes[j], string("onlynai.am.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			calc.AddCalculation(detBCf, calcDists[i], calcTypes[j], string("onlynai.cf.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			
			calc.AddCalculation(detUNPu, calcDists[i], calcTypes[j], string("uppdown.pu.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			calc.AddCalculation(detUNAm, calcDists[i], calcTypes[j], string("uppdown.am.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			calc.AddCalculation(detUNCf, calcDists[i], calcTypes[j], string("uppdown.cf.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			
			calc.AddCalculation(det10Pu, calcDists[i], calcTypes[j], string("pvc10cm.pu.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			calc.AddCalculation(det10Am, calcDists[i], calcTypes[j], string("pvc10cm.am.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
			calc.AddCalculation(det10Cf, calcDists[i], calcTypes[j], string("pvc10cm.cf.") + dist_string + calcNames[j] + string("."), fpPerHour, beta, uncertLoops);
		}
	}
	calc.RunCalculations();
}

void CalcActTest() {
	double time = 18.213217874344792;
	double beta = 0.05;
	double alpha = time / 3600.0;
	
	std::string bkgString("pvc10cm_bkg");
	std::string detString("PVC 10cm NaI");
	std::string setupString("pvc10cm");
	vector<string> distsPu{"2.5", "5", "10"};
	vector<double> calcDists{2.5, 5.0, 10.0, 20.0};
	vector<string> dists{"5", "10", "15"};
	
	cout << detString << " det. PuO-prep" << endl;
	cout << "------------------" << endl;
	cout << "Alpha:" << alpha << endl;
	cout << "Beta:" << beta << endl;
	Detector detPu(bkgString, setupString + "_avst_@m_pu", distsPu, setupString + "_@_5m_pu", 4.74, 0.1, false);
	ArrayXd x;
	ArrayXd y;
	detPu.SetDistance(20.0);
	detPu.SetIntegrationTime(time);
	detPu.CalcActivityLimits(alpha, beta, Detector::CalcType::MEAN, x, y);
	plotData("PythonPlot", x, y);
}

void plotData(string plotFunc, const ArrayXd &x, const ArrayXd &y) {
	vector<double> targetX;
	vector<double> targetY;
	
	for (int i = 0; i < x.size(); i++) {
		targetX.push_back(x[i]);
		targetY.push_back(y[i]);
	}
	plotData(plotFunc, targetX, targetY);
}

void plotData(string plotFunc, vector<double> xvec, vector<double> yvec) {
	
	PyObject *pName, *pModule, *pFunc;
	PyObject *pArgTuple, *pValue, *pXVec, *pYVec;
	
	int i;
	Py_Initialize();
	
	PyObject *sys = PyImport_ImportModule("sys");
	PyObject *path = PyObject_GetAttrString(sys, "path");
	PyList_Append(path, PyString_FromString("."));
	PyList_Append(path, PyString_FromString("/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/"));
	PyList_Append(path, PyString_FromString("/usr/local/texlive/2014//bin/x86_64-darwin/"));
	
	pName = PyString_FromString("PythonPlot");   //Get the name of the module
	pModule = PyImport_Import(pName);     //Get the module
	
	Py_DECREF(pName);
	
	if (pModule != NULL) {
		pFunc = PyObject_GetAttrString(pModule, plotFunc.c_str());   //Get the function by its name
		/* pFunc is a new reference */
		
		if (pFunc && PyCallable_Check(pFunc)) {
			
			//Set up a tuple that will contain the function arguments. In this case, the
			//function requires two tuples, so we set up a tuple of size 2.
			pArgTuple = PyTuple_New(2);
			
			//Transfer the C++ vector to a python tuple
			pXVec = PyTuple_New(xvec.size());
			for (i = 0; i < xvec.size(); ++i) {
				pValue = PyFloat_FromDouble(xvec[i]);
				if (!pValue) {
					Py_DECREF(pXVec);
					Py_DECREF(pModule);
					fprintf(stderr, "Cannot convert array value\n");
					return;
				}
				PyTuple_SetItem(pXVec, i, pValue);
			}
			
			//Transfer the other C++ vector to a python tuple
			pYVec = PyTuple_New(yvec.size());
			for (i = 0; i < yvec.size(); ++i) {
				pValue = PyFloat_FromDouble(yvec[i]);
				if (!pValue) {
					Py_DECREF(pYVec);
					Py_DECREF(pModule);
					fprintf(stderr, "Cannot convert array value\n");
					return;
				}
				PyTuple_SetItem(pYVec, i, pValue); //
			}
			
			//Set the argument tuple to contain the two input tuples
			PyTuple_SetItem(pArgTuple, 0, pXVec);
			PyTuple_SetItem(pArgTuple, 1, pYVec);
			
			//Call the python function
			pValue = PyObject_CallObject(pFunc, pArgTuple);
			
			Py_DECREF(pArgTuple);
			Py_DECREF(pXVec);
			Py_DECREF(pYVec);
			
			if (pValue != NULL) {
				printf("Result of call: %ld\n", PyInt_AsLong(pValue));
				Py_DECREF(pValue);
			}
			
			//Some error catching
			else {
				Py_DECREF(pFunc);
				Py_DECREF(pModule);
				PyErr_Print();
				fprintf(stderr,"Call failed\n");
				return;
			}
		}
		else {
			if (PyErr_Occurred())
				PyErr_Print();
			fprintf(stderr, "Cannot find function.\n");
		}
		Py_XDECREF(pFunc);
		Py_DECREF(pModule);
	}
	else {
		PyErr_Print();
		fprintf(stderr, "Failed to load.\n");
		return;
	}
	Py_Finalize();
	return;
}