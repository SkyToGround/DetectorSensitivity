//
//  Detector.cpp

#include "Detector.h"
#include "Functors.h"

Detector::Detector(BkgResponse bkg, DistResponse distResp, AngularResponse angResp, double edge_limit, unsigned int mean_iters, unsigned int sim_iters) : bkg(bkg), distResp(distResp), angResp(angResp), distance(2.0), integrationTime(1.0), velocity(8.333), edge_limit(edge_limit), mean_iters(mean_iters), sim_iters(sim_iters) {
	CalcStartStopLM();
}

Detector::Detector() : bkg(), distResp(), angResp() {
	
}

void Detector::RandomizeParameters() {
	double newBkg = bkg.GetRandomizedCPS();
	distResp.Randomize(newBkg);
	angResp.Randomize(newBkg);
}

std::pair<double, double> Detector::GetIntegrationTimes(int m, double F) const {
	return std::pair<double, double>(integrationTime * (double(m) + F), integrationTime * (double(m) + F + 1));
}

std::vector<std::pair<double,double>> Detector::GetIntTimes(CalcType tp) const {
	std::vector<std::pair<double, double>> retVec;
	double F = 0;
	if (CalcType::LIST_MODE == tp) {
		retVec.push_back(std::pair<double, double>(startTime, stopTime));
		return retVec;
	} else if (CalcType::BEST == tp) {
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

double Detector::dist_f(const double t) const {
	return sqrt((velocity * t) * (velocity * t) + distance * distance);
}

ArrayXd Detector::dist_f(const ArrayXd &t) const {
	return sqrt((velocity * t) * (velocity * t) + distance * distance);
}

double Detector::ang_f(const double t) const {
	return asin(distance / dist_f(t));
}

ArrayXd Detector::ang_f(const ArrayXd &t) const {
	return asin(distance / dist_f(t));
}

double Detector::S(const double t) const {
	double actualDist = dist_f(t);
	double distPart = distResp(actualDist);
	double angPart = angResp(ang_f(t));
	return distPart * angPart;
}

ArrayXd Detector::S(const ArrayXd &t) const {
	return distResp(dist_f(t)) * angResp(ang_f(t));
}

double Detector::S_best(const double i_time) {
	return Int_S(0, i_time / 2.0) * 2.0;
}

double Detector::S_worst(const double i_time) {
	return Int_S(0, i_time);
}

ArrayXd Detector::S_Int(const std::vector<std::pair<double, double> > i_time) const {
	ArrayXd retArr = ArrayXd::Zero(i_time.size());
	for (int i = 0; i < i_time.size(); i++) {
		retArr[i] = Int_S(i_time[i].first, i_time[i].second);
	}
	return retArr;
}

ArrayXd Detector::S_mean(const std::vector<std::pair<double, double> > i_time) const {
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

double Detector::Int_S(const double start, const double stop) const {
	const int parts = 256;
	ArrayXd times = ArrayXd::LinSpaced(parts, start, stop);
	ArrayXd values = S(times);
	double stepSize = abs(stop - start) / (parts - 1);
	return (((values.segment(0, parts - 1) + values.segment(1,  parts - 1)) / 2.0) * stepSize).sum();
}

ArrayXd Detector::S_m(const ArrayXd &t, const double i_time) const {
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

unsigned int Detector::CriticalLimitFPH(const double fph, const CalcType tp) const {
	if (CalcType::LIST_MODE == tp) {
		return CriticalLimitLM_FPH(fph);
	}
	return CriticalLimit((fph * integrationTime) / 3600.0);
}

unsigned int Detector::CriticalLimit(const double alpha) const {
	double B = simBkg * integrationTime;
	unsigned int i = 0;
	double res = 1.0;
	while (res >= alpha) {
		if (i > boost::math::max_factorial<double>::value or i > 100) {
			goto no_factorial;
		}
		res = res - ((pow(B, i))*exp(-B)) / boost::math::factorial<double>((unsigned int)i);
		i++;
	}
	return i; //Must not be i - 1 as we are integrating to C_L - 1 according to the equation
no_factorial:
	while (res >= alpha) {
		res -= exp(-B + i*log(B) - log((1.0 + 1.0/(12.0*i) + 1.0/(288.0 * i * i))*sqrt(2*pi*i))- i*log(i/e));
		i++;
	}
	return i; //Must not be i - 1 as we are integrating to C_L - 1 according to the equation
}

unsigned int Detector::CriticalLimitLM_FPH(const double fph) const {
	double alpha = 1.0 - exp(-fph);
	return CriticalLimitLM(alpha);
}

unsigned int Detector::CriticalLimitLM(const double alpha) const {
	int i = int(simBkg * integrationTime + 0.5);
	while (1.0 - F_N(i) > alpha) {
		i++;
	}
	return i;
}

void Detector::SetVelocity(double velocity) {
	Detector::velocity = velocity;
}

ArrayXd Detector::CalcSignal(double useAct, CalcType tp) {
	if (CalcType::LIST_MODE == tp) {
		double retMeanMax;
		SimMeasurements(useAct, 0, sim_iters, retMeanMax);
		ArrayXd retArr(1);
		retArr[0] = retMeanMax;
		return retArr;
	}
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
double Detector::CalcTruePositiveProbFPH(const double fph, const double testAct, CalcType tp) const {
	
	ArrayXd sig;
	std::vector<std::pair<double, double>> intTimes = GetIntTimes(tp);
	if (tp == CalcType::MEAN) {
		sig = S_mean(intTimes);
	} else if (tp == CalcType::WORST) {
		sig = S_Int(intTimes);
	} else {
		sig = S_Int(intTimes);
	}
	
	int critical_limit = CriticalLimitFPH(fph, tp);
	if (CalcType::LIST_MODE == tp) {
		double prob = SimMeasurementsOpti(testAct, critical_limit, sim_iters);
		return 1.0 - prob;
	}
	
	double B = simBkg * integrationTime;
	ArrayXd total = B + sig * testAct;
	ArrayXd res = ArrayXd::Ones(total.size());
	for (int y = 0; y < total.size(); y++) {
		for (int j = 0; j <= critical_limit - 1; j++) {
			if (j > boost::math::max_factorial<double>::value or j > 100) {
				goto no_factorial2;
			}
			res[y] -= ((pow(total[y], j))*exp(-total[y])) / boost::math::factorial<double>((unsigned int)j);
		}
		if (false) {
no_factorial2:
			res[y] = 1.0;
			for (int k = 1; k <= critical_limit - 1; k++) {
				res[y] -= exp(-total[y] + k*log(total[y]) - log((1.0 + 1.0/(12.0*k) + 1.0/(288.0*k*k))*sqrt(2*pi*k))- k*log(k/e));
			}
		}
	}
	return 1.0 - (1.0 -res).prod();
}

double Detector::CalcActivityFPH(double fph, double beta, Detector::CalcType tp) {
	if (CalcType::LIST_MODE == tp) {
		return CalcActivity(1.0 - exp(-fph), beta, tp);
	}
	return CalcActivity((fph * integrationTime) / 3600.0, beta, tp);
}

double Detector::GetSimBkg() const {
	return simBkg;
}

double Detector::CalcActivity(double alpha, double beta, CalcType tp) {
	std::vector<std::pair<double, double>> intTimes = GetIntTimes(tp);
	ArrayXd sig;
	if (CalcType::LIST_MODE == tp) {
		return CalcActivityLM(alpha, beta);
	} else if (tp == CalcType::MEAN) {
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
	
	return p[0];
}

double Detector::CalcActivityLM(double alpha, double beta) {
	double initialGuess = CalcActivity(alpha, beta, CalcType::BEST);
	unsigned int critical_limit = CriticalLimitLM(alpha);
	
	FindActivityFunctorLM functor(this, critical_limit, beta, sim_iters);
	
	VectorXd p(1);
	VectorXd res(1);
	
	//We want to find a starting value that is close to our solution
	//Fix me: are these searches really needed? (they probably are)
	double high = initialGuess;
	double low = 0.0;
	double adder = initialGuess * 0.1;
	
	double diffLimit = 0.0001;
	
	do {
		high = high + adder;
		p[0] = high;
		functor(p, res);
		adder *= 10.0;
	} while (res[0] > 0.0);
	
	
	double testPoint;

	do {
		testPoint = low + (high - low) / 2.0;
		p[0] = testPoint;
		functor(p, res);
		if (res[0] < 0) {
			high = testPoint;
		} else {
			low = testPoint;
		}
	} while (high - low > diffLimit);
	PRINT_T(std::string("Detector::CalcActivityLM(): Sim. finished with C_L = ") + lexical_cast<std::string>(critical_limit) + std::string("."));
	return p[0];
}

void Detector::SetSimBkg(double newSimBkg) {
	Detector::simBkg = newSimBkg;
}

double Detector::p_mu(const unsigned int n, const double mu) const {
	return (exp(-mu) * pow(mu, n)) / boost::math::factorial<double>(n);
}

double Detector::P_mu(const unsigned int n, const double mu) const {
	double tempRes = 0;
	for (int i = 0; i <= n; i++) {
		tempRes += pow(mu, i) / boost::math::factorial<double>(i);
	}
	return exp(-mu) * tempRes;
}

double Detector::p_mu_alt(const unsigned int n, const double mu) const {
	return (exp(-mu) * pow(mu, n)) / boost::math::factorial<double>(n);
}

double Detector::P_mu_alt(const unsigned int n, const double mu) const {
	double tempRes = 0;
	for (int i = 1; i <= n; i++) {
		tempRes += pow(mu, i) / boost::math::factorial<double>(i);
	}
	return tempRes;
}

double Detector::F_N(const unsigned int n) const {
	const double totalTime = 3600.0;
	double sup_part;
	if (n > 100) {
		sup_part = (1.0 - (simBkg * integrationTime) / (n + 1.0)) * simBkg * (totalTime - integrationTime) * p_mu_alt(n, simBkg * integrationTime);
		return  P_mu_alt(n, simBkg * integrationTime) * exp(-sup_part);
	}
	sup_part = (1.0 - (simBkg * integrationTime) / (n + 1.0)) * simBkg * (totalTime - integrationTime) * p_mu(n, simBkg * integrationTime);
	return  P_mu(n, simBkg * integrationTime) * exp(-sup_part);
}

double Detector::SimMeasurements(double actFac, unsigned int critical_limit, unsigned int iterations, double &meanMax) const {
	unsigned int truePositiveProb = 0;
	double maxRate = S(0.0) * actFac + simBkg;
	
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
	boost::random::exponential_distribution<> expDist(maxRate);
	boost::random::uniform_real_distribution<> rejectDist(0, maxRate);
	
	double cTime, pTime;
	
	//Fix me: do I really need cCount? Maybe queue::size() is fast enough?
	//This code should be profiled, boost::circular_buffer might be faster
	unsigned int cCount, maxValue = 0;
	std::queue<double> eventQueue;
	unsigned int sumMax = 0; //Used to calculate the mean max
	double cStopTime;
	for (int i = 0; i < iterations; i++) {
		maxValue = 0;
		if (startTime>-integrationTime) {
			cTime = -integrationTime + expDist(gen);
			cStopTime = integrationTime;
		} else {
			cTime = startTime + expDist(gen);
			cStopTime = stopTime;
		}
		eventQueue = std::queue<double>(); //Clear the queue
		cCount = 0;
		do {
			if (rejectDist(gen) <= S(cTime) * actFac + simBkg) {
				eventQueue.push(cTime);
				cCount++;
				pTime = cTime - integrationTime;
				while (eventQueue.front() < pTime) {
					cCount--;
					eventQueue.pop();
				}
				sumMax += cCount;
				if (cCount > maxValue) {
					maxValue = cCount;
				}
			}
			cTime += expDist(gen);
		} while (cTime < cStopTime);
		if (maxValue >= critical_limit) {
			truePositiveProb++;
		}
	}
	meanMax = double(sumMax) / double(iterations);
	return 1.0 - double(truePositiveProb) / double(iterations);
}

double Detector::SimMeasurementsOpti(double actFac, unsigned int critical_limit, unsigned int iterations) const {
	unsigned int truePositiveProb = 0;
	double maxRate = S(0.0) * actFac + simBkg;
	
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
	boost::random::exponential_distribution<> expDist(maxRate);
	boost::random::uniform_real_distribution<> rejectDist(0, maxRate);
	
	double cTime, pTime;
	
	boost::circular_buffer<double> eventQueue(critical_limit);
	
	double cStopTime;
	
	for (int i = 0; i < iterations; i++) {
		if (startTime>-integrationTime) {
			cTime = -integrationTime + expDist(gen);
			cStopTime = integrationTime;
		} else {
			cTime = startTime + expDist(gen);
			cStopTime = stopTime;
		}
		eventQueue.clear(); //Clear the queue
		do {
			if (rejectDist(gen) <= S(cTime) * actFac + simBkg) {
				eventQueue.push_back(cTime);
				pTime = cTime - integrationTime;
				while (eventQueue.front() < pTime) {
					eventQueue.pop_front();
				}
				if (eventQueue.size() >= critical_limit) {
					truePositiveProb++;
					break;
				}
			}
			cTime += expDist(gen);
		} while (cTime < cStopTime);
	}
	return 1.0 - double(truePositiveProb) / double(iterations);
}

void Detector::CalcStartStopLM() {
	double tgtVal = S(0.0) * edge_limit;
	const double maxTimeStep = 0.000001;
	
	//Find the start time
	double lowerTime = -2.0;
	double upperTime = 0.0;
	double middleTime, middleVal, lowerVal, upperVal;
	while (upperTime - lowerTime > maxTimeStep) {
		lowerVal = S(lowerTime);
		middleTime = lowerTime + (upperTime - lowerTime) / 2.0;
		middleVal = S(middleTime);
		if	(lowerVal > tgtVal) {
			lowerTime *= 2.0;
		}
		if (middleVal < tgtVal) {
			lowerTime = middleTime;
		} else if (middleVal > tgtVal) {
			upperTime = middleTime;
		}
	}
	startTime = lowerTime;
	
	//Find the stop time
	lowerTime = 0.0;
	upperTime = 2.0;
	while (upperTime - lowerTime > maxTimeStep) {
		upperVal = S(upperTime);
		middleTime = lowerTime + (upperTime - lowerTime) / 2.0;
		middleVal = S(middleTime);
		if (upperVal > tgtVal) {
			upperTime *= 2.0;
		}
		
		if (middleVal < tgtVal) {
			upperTime = middleTime;
		} else if (middleVal > tgtVal) {
			lowerTime = middleTime;
		}
	}
	stopTime = upperTime;
}
