//
//  AngularResponse.h

#ifndef __NeutronDetectorSim__AngularResponse__
#define __NeutronDetectorSim__AngularResponse__

#include <string>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <eigen3/unsupported/Eigen/LevenbergMarquardt>
#include <random>
#include <chrono>
#include <cmath>
#include "Extrap1d.h"
#include <exception>

using namespace std;
using namespace boost;

const double pi = 3.141592653589793;

class AngularResponse {
public:
	AngularResponse(const std::vector<double> pulses, const std::vector<double> livetime, const std::vector<double> angle, const double bkg_cps);
	AngularResponse operator=(const AngularResponse &setObj);
	AngularResponse();
	double operator()(const double &angle) const;
	Eigen::ArrayXd operator()(const Eigen::ArrayXd &angle) const;
	~AngularResponse();
	void Randomize(double newBkg);
private:
	bool noAngResp;
	
	std::mt19937 rand;
	
	void CreateResponseFunc();
	
	vector<double> angles;
	vector<double> measCounts;
	vector<double> measTime;
	
	Eigen::ArrayXd ang;
	Eigen::ArrayXd cps;
	boost::shared_ptr<Extrap1d> ext;
};

struct squareLawFunctor : Eigen::DenseFunctor<double> {
	Eigen::ArrayXd xVals;
	Eigen::ArrayXd yVals;
	squareLawFunctor(Eigen::ArrayXd x, Eigen::ArrayXd y) : Eigen::DenseFunctor<double> (2, 3), xVals(x), yVals(y) {
	}
	
	int operator()(const Eigen::VectorXd &p, Eigen::VectorXd &diff) const {
		diff = yVals - p[0] / (4.0 * pi * xVals.pow(p[1]));
		return 0;
	}
};

class BkgResponse {
public:
	BkgResponse(const double pulses, const double livetime);
	BkgResponse();
	double GetCPS();
	double GetRandomizedCPS();
private:
	double pulses;
	double livetime;
};

class DistResponse {
public:
	DistResponse(const std::vector<double> pulses, const std::vector<double> livetime, const std::vector<double> dist, double bkg_cps, double activity, double activity_uncertainty, bool curve_fit = false); //curve_fit should probably default to false: fix me!
	DistResponse();
	DistResponse operator=(const DistResponse &setDist);
	double operator()(const double &dist) const;
	Eigen::ArrayXd operator()(const Eigen::ArrayXd &dist) const;
	void Randomize(double newBkg);
private:
	void FitData();
	
	Eigen::ArrayXd measCounts;
	Eigen::ArrayXd measTime;
	
	Eigen::ArrayXd cpsData;
	Eigen::ArrayXd distData;
	double activity;
	double activityUncertainty;
	double p1, p2;
	bool curve_fit;
	std::mt19937 rand;
};

#endif /* defined(__NeutronDetectorSim__AngularResponse__) */
