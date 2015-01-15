//
//  AngularResponse.h

#ifndef __NeutronDetectorSim__AngularResponse__
#define __NeutronDetectorSim__AngularResponse__

#include <iostream>
#include "Measurement.h"
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

using namespace std;
using namespace boost;

const double pi = 3.141592653589793;

class AngularResponse {
public:
	AngularResponse(std::string baseName, boost::shared_ptr<BaseMeasurement> bkg, bool He3 = false);
	AngularResponse(const vector<boost::shared_ptr<BaseMeasurement>> angMeasurements, const vector<double> angles, const boost::shared_ptr<BaseMeasurement> bkg);
	AngularResponse();
	AngularResponse operator=(const AngularResponse &setObj);
	double operator()(const double &angle);
	Eigen::ArrayXd operator()(const Eigen::ArrayXd &angle);
	~AngularResponse();
	void Randomize(double newBkg);
private:
	
	mt19937 rand;
	
	void CreateResponseFunc();
	
	const string dataLoc = "/Users/jonas/Documents/Forskarstuderande/Neutronartikel/Spektra/";
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

class DistResponse {
public:
	DistResponse(std::string baseName, vector<string> dist, double activity, double activityUncertainty, boost::shared_ptr<BaseMeasurement> bkg, bool He3, bool curve_fit = false);
	DistResponse(const vector<boost::shared_ptr<BaseMeasurement>> distMeas, vector<double> dist, double activity, double activityUncertainty, const boost::shared_ptr<BaseMeasurement> bkg, bool curve_fit = false);
	DistResponse();
	DistResponse operator=(const DistResponse &setDist);
	double operator()(const double &dist);
	Eigen::ArrayXd operator()(const Eigen::ArrayXd &dist);
	void Randomize(double newBkg);
private:
	void FitData();
	const string dataLoc = "/Users/jonas/Documents/Forskarstuderande/Neutronartikel/Spektra/";
	
	Eigen::ArrayXd measCounts;
	Eigen::ArrayXd measTime;
	
	Eigen::ArrayXd cpsData;
	Eigen::ArrayXd distData;
	double activity;
	double activityUncertainty;
	double p1, p2;
	bool curve_fit;
	mt19937 rand;
};

#endif /* defined(__NeutronDetectorSim__AngularResponse__) */
