//
//  Extrap1d.cpp

#include "Extrap1d.h"


Extrap1d::Extrap1d(Eigen::ArrayXd x, Eigen::ArrayXd y) {
	Extrap1d::x = x;
	Extrap1d::y = y;
	if (x.size() != y.size()) {
		std::cout << "Extrap1d::Extrap1d(): x and y does not have the same dimensions." << std::endl;
		return;
	}
	accelerator = gsl_interp_accel_alloc();
	interpolator = gsl_spline_alloc(gsl_interp_linear, x.size());
	gsl_spline_init(interpolator, x.data(), y.data(), x.size());
	min = x[0];
	max = x[x.size() - 1];
}

Extrap1d::~Extrap1d() {
	gsl_spline_free(interpolator);
	gsl_interp_accel_free(accelerator);
}

double Extrap1d::operator()(const double &x_in) {
	if (x_in < min) {
		return y.coeff(0);
	}
	if (x_in > max) {
		return y.coeff(y.size() - 1);
	}
	return gsl_spline_eval(interpolator, x_in, accelerator);
}

Eigen::ArrayXd Extrap1d::operator()(const Eigen::ArrayXd &x) {
	Eigen::ArrayXd ret(x.size());
	for (int i = 0; i < x.size(); i++) {
		ret.coeffRef(i) = (*this)(x.coeff(i));
	}
	return ret;
}