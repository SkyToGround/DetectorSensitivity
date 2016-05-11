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
#include <boost/random.hpp>
#include <boost/circular_buffer.hpp>
#include <queue>
#include "ConsolePrint.hpp"

using namespace std;
using namespace Eigen;

const double e = 2.718281828459045;

ArrayXd factorial(const ArrayXi f);
ArrayXd pow(const double base, const ArrayXd exponent);

class Detector {
public:
	enum class CalcType {BEST, MEAN, WORST, LIST_MODE};
	Detector(BkgResponse bkg, DistResponse distResp, AngularResponse angResp, double edge_limit, unsigned int mean_iters, unsigned int sim_iters, bool mt);
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
	ArrayXd S_Int(const std::vector<std::pair<double,double>> i_time) const;
	
	/*! Calculate the mean signal for a range of values of m. The integration times provided as a parameter are not used, only the number of integration times are taken into acount when deciding which values of m to use. This means that if the detector is not symmetrical, this function will provide bad values if the number of integration times is low.
	 \param i_time Integration times. Not actually used.
	 \return The mean signal for a range of values of m.
	 */
	ArrayXd S_mean(const std::vector<std::pair<double,double>> i_time) const;
	void SetVelocity(double velocity);
	double GetVelocity() const {return velocity;};
	
	/*! Find the critical limit in number of pulses based on a given acceptable false positive rate. This function calls Detector::CriticalLimit() or Detector::CriticalLimitLM_FPH() depending on the input arguments.
	 \param fph Number of acceptable false positives per hour. Does not have to be an integer but has to be a number greater than 0.
	 \param tp The calculation type for which the critical limit should be calculated.
	 \return The critical limit. Note that the actual false positive rate might be lower than the one given as a parameter.
	 */
	unsigned int CriticalLimitFPH(const double fph, const CalcType tp) const;
	
	/*! Version of Detector::CriticalLimitFPH() for list mode calculations. Calls Detector::CrticalLimitLM().
	 \param fph Number of acceptable false positives per hour. Does not have to be an integer but has to be a number greater than 0.
	 \return The critical limit. Note that the actual false positive rate might be lower than the one given as a parameter.
	 */
	unsigned int CriticalLimitLM_FPH(const double fph) const;
	
	/*! Find the critical limit in number of pulses based on a given false positive probability (alpha). The critcial limit is calculated using the poisson distribution for small values of the critical limit and the normal distribution for large values of the critical limit. Note that if the critical limit for list mode measurements are to be calculated, use Detector::CriticalLimitLM() instead.
	 \param alpha The acceptable false positive probability.
	 \return The critical limit. Note that the actual false positive probability might be lower than the one given as a parameter though it will never be higher.
	 */
	unsigned int CriticalLimit(const double alpha) const;
	
	/*! Version of Detector::CriticalLimit() for use with list mode calculations.
	 \param alpha The acceptable false positive probability.
	 \return The critical limit. Note that the actual false positive probability might be lower than the one given as a parameter though it will never be higher.
	 */
	unsigned int CriticalLimitLM(const double alpha) const;
	
	/*! Returns the integration periods used in the calculation. If the input paramter is CalcType::MEAN, the function will return 
	 
	 */
	std::vector<std::pair<double,double>> GetIntTimes(CalcType tp) const;
	
	/*! Calculate the minimum detectable activity using the current detector settings and probabilities.
	 This member function will calculate the minimum activity required to detect a source using the provided probabilities of false positives and true positives. Note that the actual false positive probability will be lower as the value provided here is used as an "acceptabel false positive" value.
	 \param alpha The acceptable false positive probability.
	 \param beta The probability of finding the source.
	 \param tp The type of integration to be used.
	 \return The minimum required activity to find the source.
	 */
	double CalcActivity(double alpha, double beta, CalcType tp);
	
	/*! Exactly the same as Detector::CalcActivity() except this function takes an acceptable false positive probability as input instead of alpha.
	 \param fph Acceptable number of false positives per hour.
	 \param beta The probability of finding the source.
	 \param tp The type of integration to be used.
	 \return The minimum required activity to find the source.
	 */
	double CalcActivityFPH(double fph, double beta, CalcType tp);
	
	/*! The mean number of pulses in the relevant integration periods for a specific set of measurement parameters. If the calculation type is CalcType::LIST_MODE, the result is the mean maximum number of pulses between the two times given by Detector::GetIntTimes().
	 \param useAct The activity of the source measured.
	 \param tp The type of integration used.
	 \return The mean number or maximum mean number of pulses expected.
	 */
	ArrayXd CalcSignal(double useAct, CalcType tp);
	
	/*! Find the actual true positive probability for a given set of measurement parameters. If CalcActivity() is used to find the minimum detectable activity, this function should yield a value very close to 1 - beta.
	 \param alpha Acceptable false positive probability.
	 \param testAct The activity factor used.
	 \param tp The type of integration used
	 \return The true positive probability.
	 */
	double CalcTruePositiveProbFPH(const double fph, const double testAct, CalcType tp) const;
	
	
	/*! Set the number of simulations to use in the loop approximating the value of beta and the mean maximum value when performing list mode calculations.
	 \param sim_iters the number of iterations.
	 
	 */
	void SetSimIters(unsigned int sim_iters) {Detector::sim_iters = sim_iters;};
	
	ArrayXd S(const ArrayXd &t) const;
	double S(const double t) const;
	ArrayXd dist_f(const ArrayXd &t) const;
	ArrayXd ang_f(const ArrayXd &t) const;
	
	void SetSimBkg(double newSimBkg);
	
	double GetSimBkg() const;
	
	void RandomizeParameters();
private:
	DistResponse distResp;
	AngularResponse angResp;
	BkgResponse bkg;
	
	double simBkg;
	
	double velocity;
	
	double distance;
	
	double edge_limit;
	double startTime, stopTime;
	
	unsigned int mean_iters;
	unsigned int sim_iters;
	bool mt;
	
	double integrationTime;
	double dist_f(const double t) const;
	double ang_f(const double t) const;
	ArrayXd S_m(const ArrayXd &t, const double i_time) const;
	double Int_S(const double start, const double stop) const;
	std::pair<double,double> GetIntegrationTimes(int m, double F) const;
	
	/*! Poisson probability distribution function.
	 \param n Index.
	 \param mu Expected value.
	 \return Probability.
	 */
	double p_mu(const unsigned int n, const double mu) const;
	double p_mu_alt(const unsigned int n, const double mu) const;
	
	/*! Poisson cumulative probability function.
	 \param n Index.
	 \param mu Expected value.
	 \return Cumulative probability.
	 */
	double P_mu(const unsigned int n, const double mu) const;
	double P_mu_alt(const unsigned int n, const double mu) const;
	
	/*! Cumulative probablity function of max value of scanning process.
	 \param n The tested value.
	 \return The probability that the max value in the scan is equal to or less than n.
	 */
	double F_N(const unsigned int n) const;
	
	/*! List mode version of Detector::CalcActivity(). Called by Detector::CalcActivity().
	 \param alpha False positive probability per hour.
	 \param beta False negative probability.
	 \return Minimum detectable activity.
	 */
	double CalcActivityLM(double alpha, double beta);
	
	/*! Calculates the false negative probability in list mode measurements using simulations.
	 \param actFac The factor which modifies the amplitude of the source response curve.
	 \param critical_limit The limit used to decide if the source was detected or not.
	 \param iterations The number of iterations used when running the simulation.
	 \return meanMax Holds the mean maximum value from the simulations run on return.
	 \return The probability of not detecting a radioactive source which is present.
	 */
	double SimMeasurements(double actFac, unsigned int critical_limit, unsigned int iterations, double &meanMax) const;
	
	/*! Calculates the false negative probability in list mode measurements using simulations. A slightly optimized version of Detector::SimMeasurements().
	 \param actFac The factor which modifies the amplitude of the source response curve.
	 \param critical_limit The limit used to decide if the source was detected or not.
	 \param iterations The number of iterations used when running the simulation.
	 \return The probability of not detecting a radioactive source which is present.
	 */
	double SimMeasurementsOpti(double actFac, unsigned int critical_limit, unsigned int iterations) const;
	
	/* Uses binary search to find the time limits used in list mode simulations.
	 */
	void CalcStartStopLM();
	friend struct FindActivityFunctorLM;
};

#endif /* defined(__NeutronDetectorSim__Detector__) */
