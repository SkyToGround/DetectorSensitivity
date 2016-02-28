//
//  CalcCoordinator.h

#ifndef __NeutronDetectorSim__CalcCoordinator__
#define __NeutronDetectorSim__CalcCoordinator__

#include <locale>
#include <string>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <vector>
#include <iostream>
#include <ctime>
#include <boost/shared_ptr.hpp>
#include "Detector.h"
#include <utility>
#include "OutputResult.hpp"

using namespace std;

struct CalcInfo {
	CalcInfo(Detector det, Detector::CalcType cType, double fpPerHour, double beta, unsigned int uncertaintyCalcLoops = 0) : det(det), cType(cType), uncertLoops(uncertaintyCalcLoops), fpPerHour(fpPerHour), beta(beta), plotCalc(false), start_time(0.0), stop_time(0.0), steps(0), fixedInt(false), intTime(1.0) {
	}
	CalcInfo() : plotCalc(false), fixedInt(false) {
	}
	void SetPlotDetails(double start_time, double stop_time, unsigned int steps) {
		plotCalc = true;
		CalcInfo::start_time = start_time;
		CalcInfo::stop_time = stop_time;
		CalcInfo::steps = steps;
	}
	void FixIntegrationTime(double intTime) {
		fixedInt = true;
		CalcInfo::intTime = intTime;
	}
	Detector det;
	double fpPerHour;
	double beta;
	int uncertLoops;
	Detector::CalcType cType;
	bool plotCalc;
	double start_time, stop_time;
	unsigned int steps;
	bool fixedInt;
	double intTime;
};

/*! Find the global minima of a activity vs integration time function. This function does this by repeatedly calling the FindLocalMinima() for higher and higher values of C_L until the local minima starts to icnrease instead of decrease.
 \param info All the relevant measurement parameters.
 \param time Holds the optimal integration time when the function returns.
 \param bestM Holds the the value of the function at the calculated time on return.
 */
void FindGlobalMinima(CalcInfo &info, double &time, double &bestM);

/*! Find the local minima of a activity vs integration time function.
 Due to the use of Poisson-statistics, the function which describes the minimum activity which can be located by the detectorsystem for a given integration time, is a form of a sawtooth function.
 This function finds a local minima in such a function.
 \param info Contains all the relevant measurement parameters.
 \param target_C_L The critical limit for which the minima should be found.
 \param iTime The shortest possible integration time. This variable will contain the time at which the minima was found on return.
 \param minM This variable will contain the value of the function at time iTime.
 */
void FindLocalMinima(CalcInfo &info, int target_C_L, double &iTime, double &minM);

class CalcThread {
public:
	CalcThread(int threadId, boost::shared_ptr<vector<CalcInfo>> calcValues, boost::shared_ptr<boost::mutex> calcValuesMutex, boost::shared_ptr<boost::mutex> outputMutex, boost::shared_ptr<OutputResult> outResults); // boost::shared_ptr<property_tree::ptree> ptResults, boost::shared_ptr<boost::mutex> ptMutex, OutputType out_type,
	/*! This function is the main loop of the thread. It simply picks a data set to perform calculations on and calls CalcThread::CalcTimeAndAct().
	 */
	void operator()();
private:
	int CalcTimeAndAct(CalcInfo &info);
	int CalcTimeVsActData(CalcInfo & info);
	boost::shared_ptr<vector<CalcInfo>> calcValues;
	boost::shared_ptr<boost::mutex> calcValuesMutex;
	boost::shared_ptr<boost::mutex> outputMutex;
	boost::shared_ptr<OutputResult> outResults;
	int threadId;
};

class CalcCoordinator {
public:
	CalcCoordinator(OutputResult::OutputType out_t, std::string resultFileName = "results.txt");
	/*! Starts as many calculation threads as there are cores in the computer and then waits for all the threads to exit.
	 */
	void RunCalculations();
	void AddCalculation(const Detector &det, Detector::CalcType cType, double fpPerHour, double beta, int uncertLoops = 0);
	void AddFixedCalculation(const Detector &det, Detector::CalcType cType, double fpPerHour, double beta, double intTime, int uncertLoops = 0);
	
	void AddPlotCalculation(const Detector &det, Detector::CalcType cType, double fpPerHour, double beta, double start_time, double stop_time, unsigned int steps);
	~CalcCoordinator();
private:
	boost::shared_ptr<boost::mutex> outputMutex;
	boost::shared_ptr<boost::mutex> calcValuesMutex;
	string outFileName;
	boost::shared_ptr<OutputResult> outResults;
	boost::shared_ptr<vector<CalcInfo>> calcInfo;
	vector<thread> calcThreads;
};

#endif /* defined(__NeutronDetectorSim__CalcCoordinator__) */
