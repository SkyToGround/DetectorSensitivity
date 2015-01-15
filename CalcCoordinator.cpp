//
//  CalcCoordinator.cpp

#include "CalcCoordinator.h"

CalcCoordinator::CalcCoordinator(OutputType out_t, std::string resultFileName) : outFileName(resultFileName), outputMutex(new boost::mutex()), calcValuesMutex(new boost::mutex()), calcInfo(new vector<CalcInfo>()), pt(new boost::property_tree::ptree()), out_t(out_t) {
	
}

CalcCoordinator::~CalcCoordinator() {
	property_tree::write_json(outFileName, *pt.get());
}


void CalcCoordinator::RunCalculations() {
	int numberOfThreads = thread::hardware_concurrency();
	
	for (int i = 0; i < numberOfThreads; i++) {
		CalcThread calcThread(i, calcInfo, calcValuesMutex, outputMutex, out_t, pt);
		calcThreads.push_back(thread(calcThread));
	}
	
	while (calcThreads.size() > 0) {
		calcThreads[calcThreads.size() - 1].join();
		calcThreads.pop_back();
	}
}

void CalcCoordinator::AddCalculation(const Detector &det, double distance, Detector::CalcType cType, string path, double fpPerHour, double beta, int uncertLoops) {
	calcInfo->push_back(CalcInfo(det, distance, path, cType, fpPerHour, beta, uncertLoops));
}


CalcThread::CalcThread(int threadId, boost::shared_ptr<vector<CalcInfo>> calcValues, boost::shared_ptr<boost::mutex> calcValuesMutex, boost::shared_ptr<boost::mutex> outputMutex, OutputType out_t, boost::shared_ptr<property_tree::ptree> ptResults) : calcValues(calcValues), calcValuesMutex(calcValuesMutex), threadId(threadId), outputMutex(outputMutex), out_t(out_t), ptResults(ptResults) {
	
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
		CalcTimeAndAct(tmpInfo);
	}
}

int CalcThread::CalcTimeAndAct(CalcInfo &info) {
	vector<pair<string,string>> result;
	info.det.SetDistance(info.calcDist);
	double time, actValue;
	vector<double> timeVec;
	vector<double> actVec;
	
	FindGlobalMinima(info, time, actValue);
	
	result.push_back(pair<string, string>("time", lexical_cast<string>(time)));
	result.push_back(pair<string, string>("activity", lexical_cast<string>(actValue)));
	info.det.SetIntegrationTime(time);
	double srcDetectProb = info.det.CalcTruePositiveProb((info.fpPerHour * time) / 3600.0, time, info.cType);
	result.push_back(pair<string, string>("srcDetProb", lexical_cast<string>(srcDetectProb)));
	result.push_back(pair<string, string>("critical_limit", lexical_cast<string>(info.det.CriticalLimit((info.fpPerHour * time) / 3600.0))));
	result.push_back(pair<string, string>("mean_signal", lexical_cast<string>(info.det.CalcSignal(info.cType))));
	
	if (info.uncertLoops > 1) {
		double tempTime, tempAct, timeSum = 0, actSum = 0;
		for (int y = 0; y < info.uncertLoops; y++) {
			info.det.RandomizeParameters();
			FindGlobalMinima(info, tempTime, tempAct);
			actSum += tempAct;
			timeSum += tempTime;
			timeVec.push_back(tempTime);
			actVec.push_back(tempAct);
		}
		double actMean = actSum / info.uncertLoops;
		double timeMean = timeSum / info.uncertLoops;
		
		double timeDiffSum = 0;
		double actDiffSum = 0;
		
		for (int i = 0; i < info.uncertLoops; i++) {
			timeDiffSum += (timeVec[i] - timeMean) * (timeVec[i] - timeMean);
			actDiffSum += (actVec[i] - actMean) * (actVec[i] - actMean);
		}
		
		double timeStdDev = sqrt(timeDiffSum / info.uncertLoops);
		double actStdDev = sqrt(actDiffSum / info.uncertLoops);
		result.push_back(pair<string, string>("timeStdDev", lexical_cast<string>(timeStdDev)));
		result.push_back(pair<string, string>("activityStdDev", lexical_cast<string>(actStdDev)));
		result.push_back(pair<string, string>("activityMean", lexical_cast<string>(actMean)));
		result.push_back(pair<string, string>("timeMean", lexical_cast<string>(timeMean)));
	}
	
	{
		boost::mutex::scoped_lock lock(*outputMutex.get());
		for (int i = 0; i < result.size(); i++) {
			if (out_t == OutputType::SCREEN) {
				cout << result[i].first << " : " << result[i].second << endl;
			} else if (out_t == OutputType::JSON_FILE) {
				ptResults->put(info.resultPath +  result[i].first, result[i].second);
			}
		}
	}
	return 0;
}

void FindGlobalMinima2(CalcInfo &info, double &time, double &bestM) {
	double timeCurrent, mCurrent, lastTime = 0.5, lastM;
	
	FindLocalMinima(info.det, info.fpPerHour, info.beta, info.cType, lastTime, lastM);
	timeCurrent = lastTime + 0.05;
	while (true) {
		FindLocalMinima(info.det, info.fpPerHour, info.beta, info.cType, timeCurrent, mCurrent);
		if (mCurrent < lastM) {
			lastTime = timeCurrent;
			lastM = mCurrent;
		} else {
			break;
		}
		timeCurrent += 0.05;
	}
	time = lastTime;
	bestM = lastM;
}

void FindGlobalMinima(CalcInfo &info, double &time, double &bestM) {
	double currentTime = 0.0, lastTime = 5.0;
	double lastM, currentM;
	int current_C_L = 2; //Fix me, could possible need to be a larger value for systems with a very low background
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
	
	info.det.SetVelocity(8.3333); //Fix me
	info.det.SetDistance(info.calcDist);
	double alpha;
	
	//Find where we have C_L = target_C_L
	int current_C_L;
	do {
		info.det.SetIntegrationTime(upper_time);
		alpha = (info.fpPerHour * upper_time) / 3600.0;
		current_C_L = info.det.CriticalLimit(alpha);
		
		if (current_C_L > target_C_L) {
			upper_time = lower_time + (upper_time - lower_time) / 2.0;
		} else if (current_C_L < target_C_L) {
			lower_time = upper_time;
			upper_time *= 2;
		}
	 } while (current_C_L != target_C_L);
	
	lower_time = upper_time; //This is now or minimum time
	
	//Find some integration time where C_L > 1
	do {
		upper_time *= 2.0; //Potentially fix me, probably a bit to radical to double the upper time
		info.det.SetIntegrationTime(upper_time);
		alpha = (info.fpPerHour * upper_time) / 3600.0;
	} while (info.det.CriticalLimit(alpha) == target_C_L);
	
	//Now do the actual minima-finding part
	double c_time;
	while (upper_time - lower_time > maxTimeDiff) {
		c_time = (upper_time - lower_time) / 2.0 + lower_time;
		info.det.SetIntegrationTime(c_time);
		alpha = (info.fpPerHour * c_time) / 3600.0;
		current_C_L = info.det.CriticalLimit(alpha);
		if (current_C_L == target_C_L) {
			lower_time = c_time;
		} else {
			upper_time = c_time;
		}
	}
	iTime = lower_time;
	alpha = (info.fpPerHour * iTime) / 3600.0;
	info.det.SetIntegrationTime(iTime);
	minM = info.det.CalcActivity(alpha, info.beta, info.cType);
}

void FindLocalMinima(Detector &det, double falsePositivePerHour, double beta, Detector::CalcType calcType, double &iTime, double &minM) {
	const double timeStep = 0.05;	//Initial time step size, might have to be modified
	
	double lowerTime = iTime;	//The shortest possible integration time.
	double upperTime = 0;		//The variable which will contain a time which is above the valley.
	det.SetIntegrationTime(lowerTime);
	double alpha = (falsePositivePerHour * lowerTime) / 3600.0;
	double lowerM = det.CalcActivity(alpha, beta, calcType);	//Calculate a starting point
	double tempM;
	
	//cout << "Starting activity:" << lowerM << endl;
	
	//Find where we rapidly shift to higher M-value
	while (true) {
		upperTime = lowerTime + timeStep;
		det.SetIntegrationTime(upperTime);
		alpha = (falsePositivePerHour * upperTime) / 3600.0;
		tempM = det.CalcActivity(alpha, beta, calcType);
		if (tempM < lowerM) {
			lowerTime = upperTime;
			lowerM = tempM;
		} else {
			break;
		}
	}
	
	double testTime;
	//Now pinpoint where this happens
	while (true ) {
		testTime = (upperTime - lowerTime) / 2.0 + lowerTime;
		det.SetIntegrationTime(testTime);
		alpha = (falsePositivePerHour * testTime) / 3600.0;
		tempM = det.CalcActivity(alpha, beta, calcType);
		if (tempM < lowerM) {
			if (abs(1 - lowerTime / testTime) < 0.005 and abs(1 - lowerM / tempM) < 0.005) {
				minM = tempM;
				iTime = testTime;
				return;
			}
			else {
				lowerM = tempM;
				lowerTime = testTime;
			}
		} else {
			upperTime = testTime;
		}
	}
}