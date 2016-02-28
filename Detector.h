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
	ArrayXd S;
	double B;
	unsigned int C_L;
	double beta;
	FindActivityFunctor(ArrayXd S, double B, unsigned int C_L, double beta) : DenseFunctor<double> (1, 1), S(S), B(B), C_L(C_L), beta(beta) {
		
	}
	
	//Should yield the same result as CalcTruePositiveProb
	int operator()(const VectorXd &p, VectorXd &res) const {
		double M = p[0];
		ArrayXd tot = S * M + B;
		res = VectorXd(1);
		ArrayXi i = ArrayXi::LinSpaced(C_L, 0, C_L - 1);
		ArrayXd tgt = ArrayXd::Zero(tot.size());
		for (int j = 0; j < tot.size(); j++) {
			if (C_L >= boost::math::max_factorial<double>::value) {
				tgt[j] = gauss(i.cast<double>(), tot[j]).sum();
			} else {
				tgt[j] = exp(-tot[j])*(pow(tot[j], i.cast<double>())/factorial(i)).sum();
			}
		}
		res[0] = tgt.prod() - beta;
		return 0;
	}
//	int operator()(const VectorXd &p, VectorXd &res) const {
//		double M = p[0];
//		ArrayXd tot = S * M + B;
//		res = VectorXd(1);
//		ArrayXi i = ArrayXi::LinSpaced(C_L, 0, C_L - 1);
//		if (C_L >= boost::math::max_factorial<double>::value) {
//			res[0] = gauss(i.cast<double>(), tot).sum() - beta;
//		} else {
//			res[0] = exp(-tot)*(pow(tot, i.cast<double>())/factorial(i)).sum() - beta;
//		}
//		return 0;
//	}
};

class Detector {
public:
	enum class CalcType {BEST, MEAN, WORST, LIST_MODE};
	Detector(BkgResponse bkg, DistResponse distResp, AngularResponse angResp, double edge_limit, unsigned int mean_iters);
	Detector();
	~Detector();
	void SetDistance(double distance);
	double GetDistance() {return distance;};
	void SetIntegrationTime(double intTime);
	double S_best(const double i_time);
	double S_worst(const double i_time);
	double S_mean(const double i_time);
	
	/*! Integrate the response function over the integration times provided as an input. Currently uses a very simple and inaccurate integration algorithm which probably should be switched for something better.
	 \param i_time A vector of integration times to use.
	 \param Returns the integration values as Eigen::ArrayXd in order to simplify calculations later.
	 */
	ArrayXd S_Int(const std::vector<std::pair<double,double>> i_time);
	
	/*! Calculate the mean signal for a range of values of m. The integration times provided as a parameter are not used, only the number of integration times are taken into acount when deciding which values of m to use. This means that if the detector is not symmetrical, this function will provide bad values if the number of integration times is low.
	 \param i_time Integration times. Not actually used.
	 \return The mean signal for a range of values of m.
	 */
	ArrayXd S_mean(const std::vector<std::pair<double,double>> i_time);
	void SetVelocity(double velocity);
	double GetVelocity() {return velocity;};
	unsigned int CriticalLimitFPH(const double fph);
	unsigned int CriticalLimit(const double alpha);
	
	std::vector<std::pair<double,double>> GetIntTimes(CalcType tp);
	
	/*! Calculate the minimum detectable activity using the current detector settings and probabilities.
	 This member function will calculate the minimum activity required to detect a source using the provided probabilities of false positives and true positives. Note that the actual false positive probability will be lower as the value provided here is used as an "acceptabel false positive" value.
	 \param alpha The acceptable false positive probability.
	 \param beta The probability of finding the source.
	 \param tp The type of integration to be used.
	 \return The minimum required activity to find the source.
	 */
	double CalcActivity(double alpha, double beta, CalcType tp);
	
	/*! The mean number of pulses in the relevant integration periods for a specific set of measurement parameters.
	 \param useAct The activity of the source measured.
	 \param tp The type of integration used.
	 \return The mean number of pulses expected.
	 */
	ArrayXd CalcSignal(double useAct, CalcType tp);
	
	/*! Find the actual true positive probability for a given set of measurement parameters. If CalcActivity() is used to find the minimum detectable activity, this function should yield a value very close to 1 - beta.
	 \param alpha Acceptable false positive probability.
	 \param testAct The activity factor used.
	 \param tp The type of integration used
	 \return The true positive probability.
	 */
	double CalcTruePositiveProb(double alpha, double testAct, CalcType tp);
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
	
	double edge_limit;
	
	unsigned int mean_iters;
	
	double integrationTime;
	double dist_f(const double t);
	double ang_f(const double t);
	ArrayXd S_m(const ArrayXd &t, const double i_time);
	double Int_S(const double start, const double stop);
	std::pair<double,double> GetIntegrationTimes(int m, double F);
};

#endif /* defined(__NeutronDetectorSim__Detector__) */
