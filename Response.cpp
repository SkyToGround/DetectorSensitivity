//
//  AngularResponse.cpp

#include "Response.h"

AngularResponse::AngularResponse(const std::vector<double> pulses, const std::vector<double> livetime, const std::vector<double> angle, const double bkg_cps) : rand(boost::chrono::system_clock::now().time_since_epoch().count()), noAngResp(false) {
	if (pulses.size() != livetime.size() or pulses.size() != angle.size()) {
		throw length_error(std::string("AngularResponse(): Not the same number of angles, pulse values or live times."));
	}
	if (livetime.size() <= 1) {
		noAngResp = true;
		return;
	} else {
		ang = Eigen::ArrayXd(angle.size()); //Fix me: why are we using both ang and angles
		cps = Eigen::ArrayXd(angle.size());
		for (int i = 0; i < angle.size(); i++) {
			measCounts.push_back(pulses[i]);
			measTime.push_back(livetime[i]);
			cps[i] = pulses[i] / livetime[i] - bkg_cps;
			ang[i] = angle[i];
			AngularResponse::angles.push_back(angle[i]);
		}
	}
	CreateResponseFunc();
}

AngularResponse::AngularResponse() {
	noAngResp = true;
}

AngularResponse::~AngularResponse() {
}

void AngularResponse::CreateResponseFunc() {
	cps = cps / cps.maxCoeff();
	ext = boost::shared_ptr<Extrap1d>(new Extrap1d(ang, cps));
}

void AngularResponse::Randomize(double newBkg) {
	if (noAngResp) {
		return;
	}
	normal_distribution<double> angDist(0, pi / 18); //Uncertainty in angle is approximately 10 degrees or pi / 18
	
	for (int k = 0; k < angles.size(); k++) {
		ang[k] = angles[k];
		
		normal_distribution<double> pulseDist(measCounts[k], sqrt(measCounts[k]));
		
		cps[k] = pulseDist(rand) / measTime[k] - newBkg;
	}
	
	//Do not randomize the first and last angle measurement
	//We need the endpoints at 0 and 90 degrees
	if (ang.size() >= 3) {
		for (int i = 1; i < ang.size() - 1; i++) {
			ang[i] = ang[i] + angDist(rand);
		}
	}
	
	CreateResponseFunc();
}

AngularResponse AngularResponse::operator=(const AngularResponse &setObj) {
	noAngResp = setObj.noAngResp;
	if (noAngResp) {
		return *this;
	}
	angles = setObj.angles;
	measCounts = setObj.measCounts;
	measTime = setObj.measTime;
	ang = setObj.ang;
	cps = setObj.cps;
	CreateResponseFunc();
	return *this;
}

double AngularResponse::operator()(const double &angle) const {
	if (noAngResp) {
		return 1.0;
	}
	return (*ext)(abs(angle));
}

Eigen::ArrayXd AngularResponse::operator()(const Eigen::ArrayXd &angle) const {
	if (noAngResp) {
		return Eigen::ArrayXd::Ones(angle.size());
	}
	return (*ext)(angle.abs());
}

DistResponse::DistResponse(const std::vector<double> pulses, const std::vector<double> livetime, const std::vector<double> dist, double bkg_cps, double activity, double activity_uncertainty, bool curve_fit) : cpsData(dist.size()), distData(dist.size()), measCounts(dist.size()), measTime(dist.size()), activity(activity), activityUncertainty(activity_uncertainty), rand(boost::chrono::system_clock::now().time_since_epoch().count()), curve_fit(curve_fit) {
	
	if (dist.size() != pulses.size() or dist.size() != livetime.size()) {
		throw length_error(std::string("DistResponse(): Not the same number of distances, pulse values or live times."));
	}
	
	if (dist.size() == 0) {
		throw length_error(std::string("DistResponse(): No input pulses."));
	}
	
	for (int i = 0; i < dist.size(); i++) {
		distData[i] = dist[i];
		measCounts[i] = pulses[i];
		measTime[i] = livetime[i];
		cpsData[i] = (pulses[i] / livetime[i] - bkg_cps) / activity;
	}
	
	FitData();
}

DistResponse::DistResponse() : rand(boost::chrono::system_clock::now().time_since_epoch().count()), p1(0.0), p2(0.0) {
	
}

void DistResponse::FitData() {
	if (cpsData.size() == 1) {
		p2 = 2.0;
		p1 = cpsData[0] * 4.0 * pi * distData[0] * distData[0];
	} else if (cpsData.size() == 2 or not curve_fit) {
		double raiseSum = 0.0;
		int nrOfSums = 0;
		for (int out = 0; out < cpsData.size(); out++) {
			for (int in = out + 1; in < cpsData.size(); in++) {
				double distRatio = distData[out] / distData[in];
				double cpsRatio = cpsData[in] / cpsData[out];
				raiseSum +=log(pow(cpsRatio, 1.0 / log(distRatio)));
				nrOfSums++;
			}
		}
		p2 = raiseSum / double(nrOfSums);
		
		double origSum = 0.0;
		nrOfSums = 0;
		for (int y = 0; y < cpsData.size(); y++) {
			//origSum += cpsData[y] * 4.0 * pow(distData[y], p2);
			origSum += cpsData[y] * 4 * pi * pow(distData[y], p2);
			nrOfSums++;
		}
		p1 = origSum / double(nrOfSums);
	} else {
		Eigen::VectorXd p0;
		p0.setConstant(2, 1.0);
		squareLawFunctor fitFunctor(distData, cpsData);
		Eigen::DenseIndex nfev;
		Eigen::LevenbergMarquardt<squareLawFunctor>::lmdif1(fitFunctor, p0, &nfev);
		p1 = p0[0];
		p2 = p0[1];
	}
}

void DistResponse::Randomize(double newBkg) {
	Eigen::ArrayXd tempDist(distData);
	
	normal_distribution<double> distDist(0, 0.1); //Uncertainty in distance is 10 cm
	normal_distribution<double> actDist(activity, activityUncertainty);
	double newAct = actDist(rand);
	for (int o = 0; o < distData.size(); o++) {
		distData[o] = distData[o] + distDist(rand);
		
		normal_distribution<double> measDist(measCounts[o], sqrt(measCounts[o]));
		cpsData[o] = (measDist(rand) / measTime[o] - newBkg) / newAct;
	}
	
	FitData();
	
	distData = tempDist;
}

DistResponse DistResponse::operator=(const DistResponse &setDist) {
	measCounts = setDist.measCounts;
	measTime = setDist.measTime;
	cpsData = setDist.cpsData;
	distData = setDist.distData;
	activityUncertainty = setDist.activityUncertainty;
	activity = setDist.activity;
	p1 = setDist.p1;
	p2 = setDist.p2;
	return *this;
}

double DistResponse::operator()(const double &dist) const {
	return p1 / (4.0 * pi * pow(dist, p2));
}

Eigen::ArrayXd DistResponse::operator()(const Eigen::ArrayXd &dist) const {
	return p1 / (4.0 * pi * dist.pow(p2));
}

BkgResponse::BkgResponse(const double pulses, const double livetime) : pulses(pulses), livetime(livetime) {
	
}

BkgResponse::BkgResponse() : pulses(0.0), livetime(1.0) {
	
}

double BkgResponse::GetCPS() {
	return pulses / livetime;
}

double BkgResponse::GetRandomizedCPS() {
	mt19937 rand(boost::chrono::system_clock::now().time_since_epoch().count());
	
	normal_distribution<double> pulseDist(pulses, sqrt(pulses));
	
	return pulseDist(rand) / livetime;
}
