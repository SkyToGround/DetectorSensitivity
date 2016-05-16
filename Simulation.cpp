//
//  Simulation.cpp
//  MobileDetectorSim
//
//  Created by Jonas Nilsson on 2016-05-11.
//  Copyright © 2016 Jonas Nilsson. All rights reserved.
//

#include "Simulation.hpp"

ListModeSimulator::ListModeSimulator(unsigned int iterations, bool use_multiple_threads) : iterations(iterations), use_multiple_threads(use_multiple_threads) {
	int numberOfThreads = thread::hardware_concurrency();
	PRINT("Starting simulation threads.");
	for (int i = 0; i < numberOfThreads; i++) {
		SimulatorThread calcThread(i, &threadIn);
		threads.push_back(boost::thread(calcThread));
	}
}

ListModeSimulator::~ListModeSimulator() {
	if (threads.size() > 0) {
		CloseThreads();
	}
}

void ListModeSimulator::CloseThreads() {
	PRINT("Closing simulation threads.");
	SimulationParams par;
	par.exit = true;
	for (int i = 0; i < threads.size(); i++) {
		threadIn.push(par);
	}
	while (threads.size() > 0) {
		threads[threads.size() - 1].join();
		threads.pop_back();
	}
}

double ListModeSimulator::PerformSimulation(SimulationParams par) {
	concurrent_queue<ReturnData> threadOut;
	par.threadOut = &threadOut;
	if (use_multiple_threads) {
		for (int j = 0; j < threads.size() - 1; j++) {
			par.iterations = iterations / threads.size();
			threadIn.push(par);
		}
		par.iterations = iterations / threads.size() +  iterations % threads.size();
		threadIn.push(par);
	} else {
		par.iterations = iterations;
		threadIn.push(par);
	}
	unsigned int finishedCtr = 0;
	unsigned int truePositives = 0;
	ReturnData retVal;
	while (finishedCtr < iterations) {
		threadOut.wait_and_pop(retVal);
		finishedCtr += retVal.iterations;
		truePositives += retVal.value;
	}
	return 1.0 - double(truePositives) / double(finishedCtr);
}

SimulatorThread::SimulatorThread(int threadId, concurrent_queue<SimulationParams> *threadIn) : threadId(threadId), threadIn(threadIn) {
	
}

void SimulatorThread::operator()() {
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{static_cast<std::uint32_t>(now)  + threadId};
	double maxRate;
	double cTime, pTime;
	double cStartTime, cStopTime;
	unsigned int truePositiveProb = 0;
	
	SimulationParams cParams;
	ReturnData ret;
	while (true) {
		threadIn->wait_and_pop(cParams);
		if (cParams.exit) {
			break;
		}
		{
			auto dist = [&cParams] (double t) {return sqrt((cParams.velocity * t) * (cParams.velocity * t) + cParams.distance * cParams.distance);};
			auto ang = [&cParams, &dist] (double t) {return asin(cParams.distance / dist(t));};
			auto S = [&cParams, &dist, &ang] (double t) {return cParams.distResp(dist(t)) * cParams.angResp(ang(t));};
			maxRate = S(0.0) * cParams.activityFactor + cParams.background;
			boost::random::exponential_distribution<> expDist(maxRate);
			boost::random::uniform_real_distribution<> rejectDist(0, maxRate);
			boost::circular_buffer<double> eventQueue(cParams.critical_limit);
			truePositiveProb = 0;
			if (cParams.startTime > -cParams.integrationTime) {
				cStartTime = -cParams.integrationTime;
				cStopTime = cParams.integrationTime;
			} else {
				cStartTime = cParams.startTime;
				cStopTime = cParams.stopTime;
			}
			
			for (int i = 0; i < cParams.iterations; i++) {
				cTime = cStartTime + expDist(gen);
				eventQueue.clear(); //Clear the queue
				do {
					if (rejectDist(gen) <= S(cTime) * cParams.activityFactor + cParams.background) {
						eventQueue.push_back(cTime);
						pTime = cTime - cParams.integrationTime;
						while (eventQueue.front() < pTime) {
							eventQueue.pop_front();
						}
						if (eventQueue.size() >= cParams.critical_limit) {
							truePositiveProb++;
							break;
						}
					}
					cTime += expDist(gen);
				} while (cTime < cStopTime);
			}
			
		}
		ret.iterations = cParams.iterations;
		ret.value = truePositiveProb;
		cParams.threadOut->push(ret);
	}
}
