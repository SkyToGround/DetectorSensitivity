//
//  Extrap1d.h

#ifndef __NeutronDetectorSim__Extrap1d__
#define __NeutronDetectorSim__Extrap1d__

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

class Extrap1d {
public:
	Extrap1d(Eigen::ArrayXd x, Eigen::ArrayXd y);
	Extrap1d(std::vector<double> x, std::vector<double> y);
	~Extrap1d();
	double operator()(const double &x);
	Eigen::ArrayXd operator()(const Eigen::ArrayXd &x);
private:
	double min, max;
	Eigen::ArrayXd x;
	Eigen::ArrayXd y;
	gsl_spline *interpolator;
	gsl_interp_accel *accelerator;
};

#endif /* defined(__NeutronDetectorSim__Extrap1d__) */
