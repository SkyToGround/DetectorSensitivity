//
//  Detector.h

#ifndef __NeutronDetectorSim__Detector__
#define __NeutronDetectorSim__Detector__

#include <iostream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/FFT>
#include <eigen3/unsupported/Eigen/src/NonLinearOptimization/dogleg.h>
#include <eigen3/unsupported/Eigen/src/NonLinearOptimization/r1updt.h>
#include <eigen3/unsupported/Eigen/src/NonLinearOptimization/r1mpyq.h>
#include <eigen3/unsupported/Eigen/src/NonLinearOptimization/fdjac1.h>
#include <eigen3/unsupported/Eigen/src/NonLinearOptimization/HybridNonLinearSolver.h>
#include <vector>
#include "Response.h"
#include <cmath>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;
using namespace Eigen;

double gauss(const double x, const double mu);
ArrayXd gauss(const ArrayXd x, const double mu);
ArrayXd factorial(const ArrayXi f);
ArrayXd pow(const double base, const ArrayXd exponent);

struct FindActivityFunctor : DenseFunctor<double> {
	double S;
	double B;
	double C_L;
	double beta;
	FindActivityFunctor(double S, double B, int C_L, double beta) : DenseFunctor<double> (1, 1), S(S), B(B), C_L(C_L), beta(beta) {
		
	}
	
	int operator()(const VectorXd &p, VectorXd &res) const {
		double M = p[0];
		double tot = S * M + B;
		res = VectorXd(1);
		ArrayXi i = ArrayXi::LinSpaced(C_L, 0, C_L - 1);
		if (C_L >= boost::math::max_factorial<long double>::value) {
			res[0] = gauss(i.cast<double>(), tot).sum() - beta;
		} else {
			res[0] = exp(-tot)*(pow(tot, i.cast<double>())/factorial(i)).sum() - beta;
		}
		return 0;
	}
};

class Detector {
public:
	enum class CalcType {BEST, MEAN, WORST};
	Detector(BkgResponse bkg, DistResponse distResp, AngularResponse angResp, double activity);
	Detector();
	~Detector();
	void SetDistance(double distance);
	double GetDistance() {return distance;};
	void SetIntegrationTime(double intTime);
	double S_best(const double i_time);
	double S_worst(const double i_time);
	double S_mean(const double i_time);
	void SetVelocity(double velocity);
	double CriticalLimit(const double alpha);
	double CalcBoundaryTime(double alpha, int k);
	
	/*! Calculate the minimum detectable activity using the current detector settings and probabilities.
	 This member function will calculate the minimum activity required to detect a source using the provided probabilities of false positives and true positives. Note that the actual false positive probability will be lower as the value provided here is used as an "acceptabel false positive" value.
	 \param alpha The acceptable false positive probability.
	 \param beta The probability of finding the source.
	 \param tp The type of integration to be used.
	 \return The minimum required activity to find the source.
	 */
	double CalcActivity(double alpha, double beta, CalcType tp = CalcType::MEAN);
	double CalcSignal(CalcType tp);
	void CalcTimeAndActivity(double dist);
	void CalcActivityLimits(double alpha, double beta, CalcType tp, ArrayXd &x, ArrayXd &y);
	double CalcTruePositiveProb(double alpha, double i_time, CalcType tp);
	ArrayXd S(const ArrayXd &t);
	double S(const double t);
	ArrayXd dist_f(const ArrayXd &t);
	ArrayXd ang_f(const ArrayXd &t);
	
	void SetSimBkg(double newSimBkg);
	
	void RandomizeParameters();
private:
	DistResponse distResp;
	AngularResponse angResp;
	BkgResponse bkg;
	
	double simBkg;
	
	double velocity;
	
	double distance;
	
	double activity;
	
	double integrationTime;
	double dist_f(const double t);
	double ang_f(const double t);
	ArrayXd S_m(const ArrayXd &t, const double i_time);
	double Int_S(const double start, const double stop);
};

#endif /* defined(__NeutronDetectorSim__Detector__) */
