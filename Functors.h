//
//  Functors.h
//

#ifndef Functors_h
#define Functors_h

#include "Detector.h"

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
};

struct FindActivityFunctorLM : DenseFunctor<double> {
	Detector *det;
	unsigned int C_L;
	double beta;
	unsigned int iterations;
	FindActivityFunctorLM(Detector *det, unsigned int C_L, double beta, unsigned int iterations) : DenseFunctor<double> (1, 1), det(det), C_L(C_L), beta(beta), iterations(iterations) {
		
	}
	
	int operator()(const VectorXd &p, VectorXd &res) const {
		double M = p[0];
		double tempVal;
		res[0] = det->SimMeasurements(M, C_L, iterations, tempVal) - beta;
		return 0;
	}
};

#endif /* Functors_h */
