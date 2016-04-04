//
//  CalcCoordinator.cpp

#include "CalcCoordinator.h"

CalcCoordinator::CalcCoordinator(OutputResult::OutputType out_t, std::string resultFileName) : outFileName(resultFileName), outputMutex(new boost::mutex()), calcValuesMutex(new boost::mutex()), calcInfo(new vector<CalcInfo>()), outResults(new OutputResult(out_t, resultFileName)) {
	
}

CalcCoordinator::~CalcCoordinator() {
	
}


void CalcCoordinator::RunCalculations() {
	int numberOfThreads = thread::hardware_concurrency();
	
	for (int i = 0; i < numberOfThreads; i++) {
		CalcThread calcThread(i, calcInfo, calcValuesMutex, outputMutex, outResults);
		calcThreads.push_back(thread(calcThread));
	}
	
	while (calcThreads.size() > 0) {
		calcThreads[calcThreads.size() - 1].join();
		calcThreads.pop_back();
	}
}

void CalcCoordinator::AddCalculation(const Detector &det, Detector::CalcType cType, double fpPerHour, double beta, int uncertLoops) {
	calcInfo->push_back(CalcInfo(det, cType, fpPerHour, beta, uncertLoops));
}

void CalcCoordinator::AddFixedCalculation(const Detector &det, Detector::CalcType cType, double fpPerHour, double beta, double intTime, int uncertLoops) {
	CalcInfo tmp(det, cType, fpPerHour, beta, uncertLoops);
	tmp.FixIntegrationTime(intTime);
	calcInfo->push_back(tmp);
}

void CalcCoordinator::AddPlotCalculation(const Detector &det, Detector::CalcType cType, double fpPerHour, double beta, double start_time, double stop_time, unsigned int steps) {
	CalcInfo tmp(det, cType, fpPerHour, beta, 0);
	tmp.SetPlotDetails(start_time, stop_time, steps);
	calcInfo->push_back(tmp);
}


CalcThread::CalcThread(int threadId, boost::shared_ptr<vector<CalcInfo>> calcValues, boost::shared_ptr<boost::mutex> calcValuesMutex, boost::shared_ptr<boost::mutex> outputMutex, boost::shared_ptr<OutputResult> outResults) : calcValues(calcValues), calcValuesMutex(calcValuesMutex), threadId(threadId), outputMutex(outputMutex), outResults(outResults) {
	
}

void CalcThread::operator()() {
	CalcInfo tmpInfo;
	while (true) {
		{
			boost::mutex::scoped_lock lock(*calcValuesMutex.get());
			if (calcValues->size() == 0) {
				break;
			}
			tmpInfo = calcValues->at(calcValues->size() - 1);
			calcValues->pop_back();
			cout << "Thread " << threadId << " takes on one calculation. " << calcValues->size() << " are left." << endl;
		}
		if (tmpInfo.plotCalc) {
			CalcTimeVsActData(tmpInfo);
		} else {
			CalcTimeAndAct(tmpInfo);
		}
	}
}

int CalcThread::CalcTimeVsActData(CalcInfo &info) {
	std::vector<std::pair<double, double>> timeActVec;
	std::vector<double> clTime;
	std::vector<unsigned int> clValue;
	
	double currentTime = 0.0;
	double currentM;
	int current_C_L = 2;
	while (currentTime < info.stop_time) {
		FindLocalMinima(info, current_C_L, currentTime, currentM);
		if (currentTime > info.start_time and currentTime < info.stop_time) {
			clTime.push_back(currentTime);
			clValue.push_back(current_C_L);
			
			timeActVec.push_back(std::pair<double, double>(currentTime, currentM));
		}
		current_C_L++;
	}
	
	ArrayXd testTimes = ArrayXd::LinSpaced(info.steps, info.start_time, info.stop_time);
	for (int u = 0; u < testTimes.size(); u++) {
		info.det.SetIntegrationTime(testTimes[u]);
		double act = info.det.CalcActivityFPH(info.fpPerHour, info.beta, info.cType);
		timeActVec.push_back(std::pair<double, double>(testTimes[u], act));
	}
	struct c_sorter
	{
		inline bool operator() (const std::pair<double, double>& p1, const std::pair<double, double>& p2)
		{
			return (p1.first < p2.first);
		}
	};
	
	std::sort(timeActVec.begin(), timeActVec.end(), c_sorter());
	
	std::string calcTypeStr;
	if (info.cType == Detector::CalcType::BEST) {
		calcTypeStr = string("BEST");
	} else if (info.cType == Detector::CalcType::MEAN) {
		calcTypeStr = string("MEAN");
	} else if (info.cType == Detector::CalcType::WORST) {
		calcTypeStr = string("WORST");
	} else if (info.cType == Detector::CalcType::LIST_MODE) {
		calcTypeStr = string("LIST_MODE");
	} else {
		calcTypeStr = string("UNKNOWN");
	}
	
	{
		boost::mutex::scoped_lock lock(*outputMutex.get());
		outResults->StartResult();
		outResults->Write("calc_type", calcTypeStr);
		outResults->Write("distance", info.det.GetDistance());
		outResults->Write("velocity", info.det.GetVelocity());
		outResults->Write("background", info.det.GetSimBkg());
		outResults->Write("cl_values", clValue);
		outResults->Write("cl_times", clTime);
		outResults->Write("times_and_acts", timeActVec);
		outResults->EndResult();
	}
	
	return 0;
}

int CalcThread::CalcTimeAndAct(CalcInfo &info) {
	vector<pair<string,string>> result;
	
	
	double time, actValue;
	vector<double> timeVec;
	vector<double> actVec;
	if (info.fixedInt) {
		time = info.intTime;
		info.det.SetIntegrationTime(time);
		actValue = info.det.CalcActivityFPH(info.fpPerHour, info.beta, info.cType);
	} else {
		FindGlobalMinima(info, time, actValue);
	}
	
	std::string calcTypeStr;
	if (info.cType == Detector::CalcType::BEST) {
		calcTypeStr = string("BEST");
	} else if (info.cType == Detector::CalcType::MEAN) {
		calcTypeStr = string("MEAN");
	} else if (info.cType == Detector::CalcType::WORST) {
		calcTypeStr = string("WORST");
	} else if (info.cType == Detector::CalcType::LIST_MODE) {
		calcTypeStr = string("LIST_MODE");
	} else {
		calcTypeStr = string("UNKNOWN");
	}
	info.det.SetIntegrationTime(time);
	double srcDetectProb = info.det.CalcTruePositiveProbFPH(info.fpPerHour, actValue, info.cType);
	int crit_limit = info.det.CriticalLimitFPH(info.fpPerHour, info.cType);
	
	ArrayXd signalArray = info.det.CalcSignal(actValue, info.cType);
	
	std::vector<std::pair<double, double>> periodsVector = info.det.GetIntTimes(info.cType);
	
	double timeStdDev = -1;
	double actStdDev = -1;
	double actMean = -1;
	double timeMean = -1;
	
	if (info.uncertLoops > 1) {
		double tempTime, tempAct, timeSum = 0, actSum = 0;
		double timeDiffSum = 0;
		double actDiffSum = 0;
		
		if (info.fixedInt) {
			for (int y = 0; y < info.uncertLoops; y++) {
				info.det.RandomizeParameters();
				tempAct = info.det.CalcActivityFPH(info.fpPerHour, info.beta, info.cType);
				actSum += tempAct;
				timeVec.push_back(tempTime);
				actVec.push_back(tempAct);
				
			}
			actMean = actSum / info.uncertLoops;
			timeMean = info.intTime;
			
			for (int i = 0; i < info.uncertLoops; i++) {
				actDiffSum += (actVec[i] - actMean) * (actVec[i] - actMean);
			}
			
			timeStdDev = 0.0;
			actStdDev = sqrt(actDiffSum / info.uncertLoops);
			
		} else {
			for (int y = 0; y < info.uncertLoops; y++) {
				info.det.RandomizeParameters();
				FindGlobalMinima(info, tempTime, tempAct);
				actSum += tempAct;
				timeSum += tempTime;
				timeVec.push_back(tempTime);
				actVec.push_back(tempAct);
			}
			actMean = actSum / info.uncertLoops;
			timeMean = timeSum / info.uncertLoops;
			
			for (int i = 0; i < info.uncertLoops; i++) {
				timeDiffSum += (timeVec[i] - timeMean) * (timeVec[i] - timeMean);
				actDiffSum += (actVec[i] - actMean) * (actVec[i] - actMean);
			}
			
			timeStdDev = sqrt(timeDiffSum / info.uncertLoops);
			actStdDev = sqrt(actDiffSum / info.uncertLoops);
		}
	}
	
	{
		boost::mutex::scoped_lock lock(*outputMutex.get());
		outResults->StartResult();
		outResults->Write("calc_type", calcTypeStr);
		outResults->Write("distance", info.det.GetDistance());
		outResults->Write("velocity", info.det.GetVelocity());
		outResults->Write("background", info.det.GetSimBkg());
		outResults->Write("fix_int_time", info.fixedInt);
		outResults->Write("int_time", time);
		outResults->Write("min_act", actValue);
		outResults->Write("true_pos_prob", srcDetectProb);
		outResults->Write("critical_limit", crit_limit);
		outResults->Write("int_periods", periodsVector);
		outResults->Write("mean_signal", signalArray);
		if (info.uncertLoops > 1) {
			outResults->Write("int_time_std_dev", timeStdDev);
			outResults->Write("int_time_mean", timeMean);
			outResults->Write("min_act_std_dev", actStdDev);
			outResults->Write("min_act_mean", actMean);
		}
		outResults->EndResult();
	}
	return 0;
}

void FindGlobalMinima(CalcInfo &info, double &time, double &bestM) {
	double currentTime = 0.0, lastTime = 5.0;
	double lastM, currentM;
	int current_C_L = 3; //Fix me: starting C_L value influences speed quite a bit. Find better starting value!
	FindLocalMinima(info, current_C_L, lastTime, lastM);
	currentTime = lastTime;
	current_C_L++;
	while (true) {
		FindLocalMinima(info, current_C_L, currentTime, currentM);
		if (currentM < lastM) {
			lastM = currentM;
			lastTime = currentTime;
			current_C_L++;
		} else {
			break;
		}
	}
	time = lastTime;
	bestM = lastM;
}

void FindLocalMinima(CalcInfo &info, int target_C_L, double &iTime, double &minM) {
	const double maxTimeDiff = 0.00001; //The maximum time step allowed when finding minimas, might need to be changed
	double upper_time = iTime + maxTimeDiff * 2; //To assist the algorithm slightly, we add the maximum time diff
	double lower_time = 0.0;
	
	//Find where we have C_L = target_C_L
	//This is done by stepping up and down in integration time until C_L = target_C_L
	int current_C_L;
	do {
		info.det.SetIntegrationTime(upper_time);
		current_C_L = info.det.CriticalLimitFPH(info.fpPerHour, info.cType);
		
		if (current_C_L > target_C_L) {
			upper_time = lower_time + (upper_time - lower_time) / 2.0;
		} else if (current_C_L < target_C_L) {
			lower_time = upper_time;
			upper_time *= 2;
		}
	 } while (current_C_L != target_C_L);
	
	lower_time = upper_time; //This is now or minimum time
	
	//Find some integration time where C_L > target_C_L
	//This is needed in order for the algorithm to do a binary search in the later stages of the algorithm
	do {
		upper_time *= 2.0; //Potentially fix me, probably a bit to radical to double the upper time
		info.det.SetIntegrationTime(upper_time);
	} while (info.det.CriticalLimitFPH(info.fpPerHour, info.cType) == target_C_L);
	
	//Now do the actual minima-finding part
	//This is done by moving the upper and lower time limits closer to each other which results in the algorithm closing in on when target_C_L becomes target_C_L + 1
	double c_time;
	while (upper_time - lower_time > maxTimeDiff) {
		c_time = (upper_time - lower_time) / 2.0 + lower_time;
		info.det.SetIntegrationTime(c_time);
		current_C_L = info.det.CriticalLimitFPH(info.fpPerHour, info.cType);
		if (current_C_L == target_C_L) {
			lower_time = c_time;
		} else {
			upper_time = c_time;
		}
	}
	//The actual integration time is the lower one as we really want C_L = target_C_L
	iTime = lower_time;
	info.det.SetIntegrationTime(iTime);
	minM = info.det.CalcActivityFPH(info.fpPerHour, info.beta, info.cType);
}
