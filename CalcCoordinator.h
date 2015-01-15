//
//  CalcCoordinator.h

#ifndef __NeutronDetectorSim__CalcCoordinator__
#define __NeutronDetectorSim__CalcCoordinator__

#include <locale>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <vector>
#include <iostream>
#include <ctime>
#include <boost/shared_ptr.hpp>
#include "Detector.h"
#include <utility>

using namespace std;

enum class OutputType {SCREEN, TEXT_FILE, JSON_FILE};

/*! Find the local minima of a activity vs integration time function.
  Due to the use of Poisson-statistics, the function which describes the minimum activity which can be located by the detectorsystem for a given integration time, is a form of a sawtooth function.
  This function finds a local minima in such a function.
 \param det The detector/source used as abase for the calculations.
 \param falsePositivePerHour Number of acceptable false positives per hour. Note that the number of actual false positives will be lower.
 \param beta The probability of locating a source of the calculated activity.
 \param calcType An enum which sets the type of integration to use. Worst, best or mean time scenario.
 \param iTime The shortest possible integration time. This variable will contain the time at which the minima was found on return.
 \param minM This variable will contain the value of the function at time iTime.
 */
void FindLocalMinima(Detector &det, double falsePositivePerHour, double beta, Detector::CalcType calcType, double &iTime, double &minM);

struct CalcInfo {
	CalcInfo(Detector det, double calcDist, string resultPath, Detector::CalcType cType, double fpPerHour, double beta, unsigned int uncertaintyCalcLoops = 0) : det(det), calcDist(calcDist), resultPath(resultPath), cType(cType), uncertLoops(uncertaintyCalcLoops), fpPerHour(fpPerHour), beta(beta) {
	}
	CalcInfo() {
	}
	CalcInfo operator=(const CalcInfo &setData) {
		det = setData.det;
		calcDist = setData.calcDist;
		cType = setData.cType;
		resultPath = setData.resultPath;
		uncertLoops = setData.uncertLoops;
		fpPerHour = setData.fpPerHour;
		beta = setData.beta;
		
		return *this;
	}
	Detector det;
	double fpPerHour;
	double beta;
	double calcDist;
	string resultPath;
	int uncertLoops;
	Detector::CalcType cType;
};

void FindGlobalMinima(CalcInfo &info, double &time, double &bestM);

void FindLocalMinima(CalcInfo &info, int target_C_L, double &iTime, double &minM);

void FindGlobalMinima2(CalcInfo &info, double &time, double &bestM);

class CalcThread {
public:
	CalcThread(int threadId, boost::shared_ptr<vector<CalcInfo>> calcValues, boost::shared_ptr<boost::mutex> calcValuesMutex, boost::shared_ptr<boost::mutex> outputMutex, OutputType out_t, boost::shared_ptr<property_tree::ptree> ptResults); // boost::shared_ptr<property_tree::ptree> ptResults, boost::shared_ptr<boost::mutex> ptMutex, OutputType out_type,
	void operator()();
private:
	int CalcTimeAndAct(CalcInfo &info);
	int CalcTimeAndAct(CalcInfo &info, double mGuess);
	boost::shared_ptr<vector<CalcInfo>> calcValues;
	boost::shared_ptr<boost::mutex> calcValuesMutex;
	boost::shared_ptr<boost::mutex> outputMutex;
	OutputType out_t;
	boost::shared_ptr<property_tree::ptree> ptResults;
	int threadId;
};

class CalcCoordinator {
public:
	CalcCoordinator(OutputType out_t = OutputType::SCREEN, std::string resultFileName = "results.txt");
	void RunCalculations();
	void AddCalculation(const Detector &det, double distance, Detector::CalcType cType, string path, double fpPerHOur, double beta, int uncertLoops = 0);
	~CalcCoordinator();
private:
	boost::shared_ptr<boost::mutex> outputMutex;
	boost::shared_ptr<boost::mutex> calcValuesMutex;
	string outFileName;
	boost::shared_ptr<boost::property_tree::ptree> pt;
	OutputType out_t;
	boost::shared_ptr<vector<CalcInfo>> calcInfo;
	vector<thread> calcThreads;
};

#endif /* defined(__NeutronDetectorSim__CalcCoordinator__) */
