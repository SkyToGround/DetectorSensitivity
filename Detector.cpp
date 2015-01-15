//
//  Detector.cpp

#include "Detector.h"

Detector::Detector(string bkgName, string distName, vector<string> dists, string angName, double activity, double activityUncertainty, bool He3) : bkg(), distResp(), angResp(), distance(2.5), integrationTime(1.0), velocity(8.333), activity(activity) {
	if (He3) {
		bkg = boost::shared_ptr<BaseMeasurement>(new He3_Measurement(dataLoc + bkgName + string(".spc")));
	} else {
		bkg = boost::shared_ptr<BaseMeasurement>(new NaI_Measurement(dataLoc + bkgName + string(".Spc")));
	}
	
	distResp = DistResponse(distName, dists, activity, activityUncertainty, bkg, He3, true);
	angResp = AngularResponse(angName, bkg, He3);
}

Detector::Detector(boost::shared_ptr<BaseMeasurement> bkg, vector<boost::shared_ptr<BaseMeasurement>> distMeasVec, vector<double> dists, vector<boost::shared_ptr<BaseMeasurement>> angMeasVec, vector<double> angles, double activity, double activityUncertainty) : bkg(bkg), distResp(), angResp(), distance(2.5), integrationTime(1.0), velocity(8.333), activity(activity) {
	
	distResp = DistResponse(distMeasVec, dists, activity, activityUncertainty, bkg);
	angResp = AngularResponse(angMeasVec, angles, bkg);
}

Detector::Detector() : bkg(), distResp(), angResp(), activity(-1.0) {
	
}

Detector::Detector(const Detector &det) : bkg(det.bkg), distResp(det.distResp), angResp(det.angResp), distance(det.distance), integrationTime(det.integrationTime), velocity(det.velocity), activity(det.activity) {
	
}

Detector Detector::operator=(const Detector &setDet) {
	bkg = setDet.bkg;
	
	distResp = setDet.distResp;
	angResp = setDet.angResp;
	
	distance = setDet.distance;
	
	integrationTime = setDet.integrationTime;
	velocity = setDet.velocity;
	activity = setDet.activity;
	//Fix me, check me (is this correct? are we missing something?)
	return *this;
}

void Detector::RandomizeParameters() {
	double newBkg = bkg->GetRandomizedCPS();
	distResp.Randomize(newBkg);
	angResp.Randomize(newBkg);
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

double Detector::CriticalLimit(const double alpha) {
	long double B = bkg->GetCPS() * integrationTime;
	unsigned long long int i = 0;
	long double res = 1.0;
	while (res >= alpha) {
		if (i > boost::math::max_factorial<long double>::value) {
			cout << "Warning, reached a max factorial input of: " << i << endl;
			goto no_factorial;
		}
		res = res - ((pow(B, i))*exp(-B)) / boost::math::factorial<long double>((unsigned int)i);
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

double Detector::CalcBoundaryTime(double alpha, int k) {
	boost::math::chi_squared dist(2*(k+1));
	return (boost::math::quantile(dist, 1 - alpha) / 2.0) / bkg->GetCPS();
}

void Detector::SetVelocity(double velocity) {
	Detector::velocity = velocity;
}

double Detector::CalcSignal(CalcType tp) {
	double sig;
	if (tp == CalcType::MEAN) {
		sig = S_mean(integrationTime);
	} else if (tp == CalcType::WORST) {
		sig = S_worst(integrationTime);
	} else {
		sig = S_best(integrationTime);
	}
	double B = bkg->GetCPS() * integrationTime;
	return B + sig * activity;
}

double Detector::CalcTruePositiveProb(double alpha, double i_time, CalcType tp) {
	double sig;
	if (tp == CalcType::MEAN) {
		sig = S_mean(i_time);
	} else if (tp == CalcType::WORST) {
		sig = S_worst(i_time);
	} else {
		sig = S_best(i_time);
	}
	
	int critical_limit = int(CriticalLimit(alpha));
	
	
	double B = bkg->GetCPS() * i_time;
	double total = B + sig * activity;
	
	long double res = 1.0;
	for (int j = 0; j <= critical_limit; j++) {
		if (j > boost::math::max_factorial<long double>::value) {
			cout << "Warning, reached a max factorial input of: " << j << endl;
			goto no_factorial2;
		}
		res -= ((pow(total, j))*exp(-total)) / boost::math::factorial<long double>((unsigned int)j);
	}
	return res;
no_factorial2:
	res = 1.0;
	for (int k = 0; k <= critical_limit; k++) {
		res -= gauss(k, total);
	}
	return res; //Fix me, double check for off by one errors!
}

double Detector::CalcActivity(double alpha, double beta, CalcType tp) {
	double sig;
	if (tp == CalcType::MEAN) {
		sig = S_mean(integrationTime);
	} else if (tp == CalcType::WORST) {
		sig = S_worst(integrationTime);
	} else {
		sig = S_best(integrationTime);
	}
	double critical_limit = CriticalLimit(alpha);
	
	FindActivityFunctor functor(sig, bkg->GetCPS() * integrationTime, critical_limit, beta);
	
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

void Detector::CalcActivityLimits(double alpha, double beta, CalcType tp, ArrayXd &x, ArrayXd &y) {
	double sig;
	if (tp == CalcType::MEAN) {
		sig = S_mean(integrationTime);
	} else if (tp == CalcType::WORST) {
		sig = S_worst(integrationTime);
	} else {
		sig = S_best(integrationTime);
	}
	double critical_limit = CriticalLimit(alpha);
	FindActivityFunctor functor(sig, bkg->GetCPS() * integrationTime, critical_limit, beta);
	
	x = ArrayXd::LinSpaced(1000, 0.01, 2500);
	y = ArrayXd::Zero(1000);
	VectorXd out(1);
	VectorXd in(1);
	for (int i = 0; i < 1000; i++) {
		in[0] = x[i];
		functor(in, out);
		y[i] = out[0];
	}
}
