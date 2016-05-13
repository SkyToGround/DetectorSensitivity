//
//  Simulation.hpp
//  MobileDetectorSim
//
//  Created by Jonas Nilsson on 2016-05-11.
//  Copyright © 2016 Jonas Nilsson. All rights reserved.
//

#ifndef Simulation_hpp
#define Simulation_hpp

#include <vector>
#include <boost/thread/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/circular_buffer.hpp>
#include "ConsolePrint.hpp"
#include "Response.h"
#include "concurrent_queue.h"

struct ReturnData {
	unsigned int value;
	unsigned int iterations;
};

struct SimulationParams {
	SimulationParams() : exit(false) {};
	double distance;
	double velocity;
	AngularResponse angResp;
	DistResponse distResp;
	double activityFactor;
	double background;
	unsigned int critical_limit;
	double startTime;
	double stopTime;
	double integrationTime;
	
	unsigned int iterations;
	bool exit;
	concurrent_queue<ReturnData> *threadOut;
};

class SimulatorThread {
public:
	SimulatorThread(int threadId, concurrent_queue<SimulationParams> *threadIn);
	void operator()();
private:
	concurrent_queue<SimulationParams> *threadIn;
	int threadId;
};

class ListModeSimulator {
public:
	ListModeSimulator(unsigned int iterations, bool use_multiple_threads);
	~ListModeSimulator();
	double PerformSimulation(SimulationParams params);
	void CloseThreads();
private:
	unsigned int iterations;
	bool use_multiple_threads;
	concurrent_queue<SimulationParams> threadIn;
	std::vector<boost::thread> threads;
};

#endif /* Simulation_hpp */
