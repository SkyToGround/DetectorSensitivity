//
//  Detector.cpp

#include "Detector.h"

Detector::Detector(BkgResponse bkg, DistResponse distResp, AngularResponse angResp, double edge_limit, unsigned int mean_iters) : bkg(bkg), distResp(distResp), angResp(angResp), distance(2.0), integrationTime(1.0), velocity(8.333), edge_limit(edge_limit), mean_iters(mean_iters) {
}

Detector::Detector() : bkg(), distResp(), angResp() {
	
}

void Detector::RandomizeParameters() {
	double newBkg = bkg.GetRandomizedCPS();
	distResp.Randomize(newBkg);
	angResp.Randomize(newBkg);
}

std::pair<double, double> Detector::GetIntegrationTimes(int m, double F) {
	return std::pair<double, double>(integrationTime * (double(m) + F), integrationTime * (double(m) + F + 1));
}

std::vector<std::pair<double,double>> Detector::GetIntTimes(CalcType tp) {
	std::vector<std::pair<double, double>> retVec;
	double F = 0;
	if (CalcType::BEST == tp) {
		F = -0.5;
	} else if (CalcType::WORST == tp) {
		F = 0.0;
	} else { //Used when calc type is MEAN and thus gives incorrect times in this case
		F = -0.25;
	}
	std::pair<double,double> cTimes;
	int m = 0, lower_m, upper_m;
	double sigLimit = S(0.0) * edge_limit;
	while (true) {
		cTimes = GetIntegrationTimes(m, F);
		if (S(cTimes.first) > sigLimit) {
			m--;
		} else {
			lower_m = m;
			break;
		}
	}
	m = 0;
	while (true) {
		cTimes = GetIntegrationTimes(m, F);
		if (S(cTimes.second) > sigLimit) {
			m++;
		} else {
			upper_m = m;
			break;
		}
	}
	for (int y = lower_m; y <= upper_m; y++) {
		retVec.push_back(GetIntegrationTimes(y, F));
	}
	return retVec;
}

Detector::~Detector() {
}

void Detector::SetDistance(double distance) {
	Detector::distance = distance;
}

void Detector::SetIntegrationTime(double intTime) {
	Detector::integrationTime = intTime;
}

double Detector::dist_f(const double t) {
	return sqrt((velocity * t) * (velocity * t) + distance * distance);
}

ArrayXd Detector::dist_f(const ArrayXd &t) {
	return sqrt((velocity * t) * (velocity * t) + distance * distance);
}

double Detector::ang_f(const double t) {
	return asin(distance / dist_f(t));
}

ArrayXd Detector::ang_f(const ArrayXd &t) {
	return asin(distance / dist_f(t));
}

double Detector::S(const double t) {
	double actualDist = dist_f(t);
	double distPart = distResp(actualDist);
	double angPart = angResp(ang_f(t));
	return distPart * angPart;
}

ArrayXd Detector::S(const ArrayXd &t) {
	return distResp(dist_f(t)) * angResp(ang_f(t));
}

double Detector::S_best(const double i_time) {
	return Int_S(0, i_time / 2.0) * 2.0;
}

double Detector::S_worst(const double i_time) {
	return Int_S(0, i_time);
}

ArrayXd Detector::S_Int(const std::vector<std::pair<double, double> > i_time) {
	ArrayXd retArr = ArrayXd::Zero(i_time.size());
	for (int i = 0; i < i_time.size(); i++) {
		retArr[i] = Int_S(i_time[i].first, i_time[i].second);
	}
	return retArr;
}

ArrayXd Detector::S_mean(const std::vector<std::pair<double, double> > i_time) {
	int low_m = 0, high_m = 0;
	int tmpCtr = int(i_time.size()) - 1;
	while (0 < tmpCtr) {
		if (tmpCtr >= 2) {
			high_m++;
			low_m--;
			tmpCtr -= 2;
		} else if (tmpCtr >= 1) {
			low_m--;
			tmpCtr--;
		}
	}
	
	for (int m = low_m; m <= high_m; m++) {
		
	}
	ArrayXd retArr = ArrayXd::Zero(i_time.size());
	for (int i = 0; i < i_time.size(); i++) {
		retArr[i] = Int_S(i_time[i].first, i_time[i].second);
	}
	return retArr;
}

double Detector::Int_S(const double start, const double stop) {
	const int parts = 256;
	ArrayXd times = ArrayXd::LinSpaced(parts, start, stop);
	ArrayXd values = S(times);
	double stepSize = abs(stop - start) / (parts - 1);
	return (((values.segment(0, parts - 1) + values.segment(1,  parts - 1)) / 2.0) * stepSize).sum();
}

ArrayXd Detector::S_m(const ArrayXd &t, const double i_time) {
	VectorXd signal(S(t));
	double delta_size = abs(t[0] - t[1]);
	int elements = int(i_time / delta_size + 0.5);
	VectorXd kernel = VectorXd::Zero(t.size());
	for (int i  = 0; i < elements; i++) {
		kernel[i] = delta_size;
	}
	
	FFT<double> fft;
	VectorXcd freq_kernel;
	VectorXcd freq_signal;
	fft.fwd(freq_signal, signal);
	fft.fwd(freq_kernel, kernel);
	VectorXcd temp = VectorXcd(freq_kernel.array() * freq_signal.array());
	VectorXd result;
	fft.inv(result, temp);
	return result.array();
}

double Detector::S_mean(const double i_time) {
	ArrayXd sig_time = ArrayXd::LinSpaced(4096 * 2, -i_time * 2 * 2, i_time * 2 * 2); //Temp buggfix
	ArrayXd sig = S_m(sig_time, i_time);
	
	int zero_pos;
	sig.maxCoeff(&zero_pos);
	double time_delta = abs(sig_time[0] - sig_time[1]);
	int time_steps = int((i_time * 0.5) / time_delta + 0.5);
	double int_S = (sig.segment(zero_pos, time_steps).sum() * time_delta) / (i_time * 0.5); //Bug här: ibland så är zero_pos nära kanten på sig, varför?
	return int_S;
}

double gauss(const double x, const double mu) {
	double sig = sqrt(mu);
	return exp(-((x - mu) * (x - mu)) / (2.0*(sig * sig))) / (sig*2.5066282746310002);
}

ArrayXd gauss(const ArrayXd x, const double mu) {
	double sig = sqrt(mu);
	return (-((x - mu) * (x - mu)) / (2.0*(sig*sig))).exp() / (sig*2.5066282746310002);
}

ArrayXd factorial(const ArrayXi f) {
	ArrayXd retArr(f.size());
	for (int i = 0; i < f.size(); i++) {
		retArr.coeffRef(i) = boost::math::factorial<double>(f[i]);
	}
	return retArr;
}

ArrayXd pow(const double base, const ArrayXd exponent) {
	ArrayXd ret(exponent.size());
	for (int i = 0; i < exponent.size(); i++) {
		ret.coeffRef(i) = pow(base, exponent[i]);
	}
	return ret;
}

unsigned int Detector::CriticalLimitFPH(const double fph) {
	return CriticalLimit((fph * integrationTime) / 3600.0);
}

unsigned int Detector::CriticalLimit(const double alpha) {
	double B = simBkg * integrationTime;
	unsigned int i = 0;
	double res = 1.0;
	while (res >= alpha) {
		if (i > boost::math::max_factorial<double>::value) {
			cout << "Warning, reached a max factorial input of: " << i << endl;
			goto no_factorial;
		}
		res = res - ((pow(B, i))*exp(-B)) / boost::math::factorial<double>((unsigned int)i);
		i++;
	}
	return i; //Must not be i - 1 as we are integrating to C_L - 1 according to the equation
no_factorial:
	i = 0;
	res = 1.0;
	while (res >= alpha) {
		res = res - gauss(i, B);
		i++;
	}
	return i; //Must not be i - 1 as we are integrating to C_L - 1 according to the equation
}

void Detector::SetVelocity(double velocity) {
	Detector::velocity = velocity;
}

ArrayXd Detector::CalcSignal(double useAct, CalcType tp) {
	ArrayXd sig;
	std::vector<std::pair<double, double>> intTimes = GetIntTimes(tp);
	if (tp == CalcType::MEAN) {
		sig = S_mean(intTimes);
	} else if (tp == CalcType::WORST) {
		sig = S_Int(intTimes);
	} else {
		sig = S_Int(intTimes);
	}
	double B = simBkg * integrationTime;
	return B + sig * useAct;
}


//Should yield the same result as FindActivityFunctor
double Detector::CalcTruePositiveProb(double alpha, double testAct, CalcType tp) {
	ArrayXd sig;
	std::vector<std::pair<double, double>> intTimes = GetIntTimes(tp);
	if (tp == CalcType::MEAN) {
		sig = S_mean(intTimes);
	} else if (tp == CalcType::WORST) {
		sig = S_Int(intTimes);
	} else {
		sig = S_Int(intTimes);
	}
	
	int critical_limit = CriticalLimit(alpha);
	
	
	double B = simBkg * integrationTime;
	ArrayXd total = B + sig * testAct;
	ArrayXd res = ArrayXd::Ones(total.size());
	for (int y = 0; y < total.size(); y++) {
		for (int j = 0; j <= critical_limit - 1; j++) {
			if (j > boost::math::max_factorial<double>::value) {
				cout << "Warning, reached a max factorial input of: " << j << endl;
				goto no_factorial2;
			}
			res[y] -= ((pow(total[y], j))*exp(-total[y])) / boost::math::factorial<double>((unsigned int)j);
		}
		if (false) {
no_factorial2:
			res[y] = 1.0;
			for (int k = 0; k <= critical_limit - 1; k++) {
				res[y] -= gauss(k, total[y]);
			}
		}
	}
	return 1.0 - (1.0 -res).prod();
}

double Detector::CalcActivity(double alpha, double beta, CalcType tp) {
	std::vector<std::pair<double, double>> intTimes = GetIntTimes(tp);
	ArrayXd sig;
	if (tp == CalcType::MEAN) {
		sig = S_mean(intTimes);
	} else if (tp == CalcType::WORST) {
		sig = S_Int(intTimes);
	} else {
		sig = S_Int(intTimes);
	}
	unsigned int critical_limit = CriticalLimit(alpha);
	
	FindActivityFunctor functor(sig, simBkg * integrationTime, critical_limit, beta);
	
	VectorXd p(1);
	VectorXd res(1);
	
	//We want to find a starting value that is close to our solution
	double low = 0.0001;
	double upp = 0.0;
	double adder = 0.0001;
	double negativeLimit = -beta / 2.0;
	double positiveLimit = (1.0 - beta) / 2.0;
	res[0] = positiveLimit; //We dont want to add to the value of low on the first loop
	//First locate where we go to negative values
	while(res[0] > negativeLimit) {
		upp = low + adder;
		p[0] = upp;
		functor(p, res);
		if (res[0] > positiveLimit) {
			low = low + adder;
		}
		adder *= 10.0;
	}
	
	double testPoint = low + (upp - low) / 2.0;
	p[0] = testPoint;
	functor(p, res);
	// Now locate where the function falls within the positive and negative limit
	while (res[0] < negativeLimit or res[0] > positiveLimit) {
		if (res[0] >= positiveLimit) {
			low = testPoint;
		} else {
			upp = testPoint;
		}
		testPoint = low + (upp - low) / 2.0;
		p[0] = testPoint;
		functor(p, res);
	}
		
	p.setConstant(1, testPoint);
	
	// do the computation
	HybridNonLinearSolver<FindActivityFunctor> solver(functor);
	solver.hybrd1(p);
	
	//calculate final value and double check that it is close to 0
	functor(p, res);
	
	if (abs(res[0]) > 0.01) {
		cout << "Failed to find 0!" << endl;
	}
	
	return p[0];
}

void Detector::SetSimBkg(double newSimBkg) {
	Detector::simBkg = newSimBkg;
}
