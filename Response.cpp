//
//  AngularResponse.cpp

#include "Response.h"

AngularResponse::AngularResponse(std::string baseName, boost::shared_ptr<BaseMeasurement> bkg, bool He3) : ang(3), cps(3), rand(chrono::system_clock::now().time_since_epoch().count()) {
	vector<string> names;
	string tmpString = baseName;
	if (baseName == string("none")) {
		ang[0] = 0.0;
		ang[1] = pi / 4.0;
		ang[2] = pi / 2.0;
		cps[0] = 1.0;
		cps[1] = 1.0;
		cps[2] = 1.0;
		CreateResponseFunc();
		return;
	}
	names.push_back(tmpString.replace(tmpString.find("@"), 1, "90grader"));
	tmpString = baseName;
	names.push_back(tmpString.replace(tmpString.find("@"), 1, "45grader"));
	tmpString = baseName;
	names.push_back(tmpString.replace(tmpString.find("@"), 1, "avst"));
	
	angles.push_back(0.0);
	angles.push_back(pi / 4.0);
	angles.push_back(pi / 2.0);
	
	for (int i = 0; i < names.size(); i++) {
		BaseMeasurement *tmpMeas = nullptr;
		if (He3) {
			tmpMeas = new He3_Measurement(dataLoc + names[i] + ".spc");
		} else {
			tmpMeas = new NaI_Measurement(dataLoc + names[i] + ".Spc");
		}
		
//		ang[i] = angles[i];
		measCounts.push_back(tmpMeas->GetTotalPulses());
		measTime.push_back(tmpMeas->GetLiveTime());
//		cps[i] = tmpMeas->GetCPS() - bkg->GetCPS();
		
		delete tmpMeas;
	}
	
	for (int j = 0; j < angles.size(); j++) {
		ang[j] = angles[j];
		cps[j] = (measCounts[j] / measTime[j]) - bkg->GetCPS();
	}
	CreateResponseFunc();
}

AngularResponse::AngularResponse(const vector<boost::shared_ptr<BaseMeasurement>> angMeasurements, const vector<double> angles, const boost::shared_ptr<BaseMeasurement> bkg) : rand(chrono::system_clock::now().time_since_epoch().count()) {
	if (angMeasurements.size() <= 1) {
		ang = Eigen::ArrayXd(2);
		cps = Eigen::ArrayXd(2);
		ang[0] = 0.0;
		ang[1] = pi / 2.0;
		cps[0] = 1.0;
		cps[1] = 1.0;
		CreateResponseFunc();
		return;
	}
	if (angMeasurements.size() != angles.size()) {
		cout << "Error: Number of angles and angle measurements are not equal!" << endl;
	}
	
	ang = Eigen::ArrayXd(angles.size());
	cps = Eigen::ArrayXd(angles.size());
	for (int i = 0; i < angles.size(); i++) {
		measCounts.push_back(angMeasurements[i]->GetTotalPulses());
		measTime.push_back(angMeasurements[i]->GetLiveTime());
		cps[i] = angMeasurements[i]->GetCPS();
		ang[i] = angles[i];
		AngularResponse::angles.push_back(angles[i]);
	}
	CreateResponseFunc();
}

AngularResponse::~AngularResponse() {
}

void AngularResponse::CreateResponseFunc() {
	cps = cps / cps.maxCoeff();
	ext = boost::shared_ptr<Extrap1d>(new Extrap1d(ang, cps));
}

void AngularResponse::Randomize(double newBkg) {
	
	normal_distribution<double> angDist(0, pi / 18); //Uncertainty in angle is approsimately 10 degrees or pi / 18
	
	for (int k = 0; k < angles.size(); k++) {
		ang[k] = angles[k];
		
		normal_distribution<double> pulseDist(measCounts[k], sqrt(measCounts[k]));
		
		cps[k] = pulseDist(rand) / measTime[k] - newBkg;
	}
	ang[1] = ang[1] + angDist(rand);
	
	CreateResponseFunc();
}

AngularResponse AngularResponse::operator=(const AngularResponse &setObj) {
	angles = setObj.angles;
	measCounts = setObj.measCounts;
	measTime = setObj.measTime;
	ang = setObj.ang;
	cps = setObj.cps;
	CreateResponseFunc();
	return *this;
}

AngularResponse::AngularResponse() : rand(chrono::system_clock::now().time_since_epoch().count()) {
	
}

double AngularResponse::operator()(const double &angle) {
	return (*ext)(abs(angle));
}

Eigen::ArrayXd AngularResponse::operator()(const Eigen::ArrayXd &angle) {
	return (*ext)(angle.abs());
}

DistResponse::DistResponse(std::string baseName, vector<string> dist, double activity, double activityUncertainty, boost::shared_ptr<BaseMeasurement> bkg, bool He3, bool curve_fit) : cpsData(dist.size()), distData(dist.size()), measCounts(dist.size()), measTime(dist.size()), activity(activity), activityUncertainty(activityUncertainty), rand(chrono::system_clock::now().time_since_epoch().count()), curve_fit(curve_fit) {
	
	if (dist.size() == 2) {
		double r1, r2;
		string c_name1 = baseName, c_name2 = baseName;
		c_name1.replace(c_name1.find("@"), 1, dist[0]);
		c_name2.replace(c_name2.find("@"), 1, dist[1]);
		NaI_Measurement m1(string("/Users/jonas/Documents/Forskarstuderande/Helsingforsresa_nr2/Helsinki/") + c_name1 + string(".Spc"));
		NaI_Measurement m2(string("/Users/jonas/Documents/Forskarstuderande/Helsingforsresa_nr2/Helsinki/") + c_name2 + string(".Spc"));
		double rad1 = lexical_cast<double>(dist[0]);
		double rad2 = lexical_cast<double>(dist[1]);
		r1 = m1.GetCPS() * 4.0 * pi * rad1 * rad1;
		r2 = m2.GetCPS() * 4.0 * pi * rad2 * rad2;
		p1 = (r1 + r2) / 2.0;
		p2 = 2.0;
		return;
	}
	
	for (int i = 0; i < dist.size(); i++) {
		BaseMeasurement *tmpMeas = nullptr;
		string c_name = baseName;
		c_name.replace(c_name.find("@"), 1, dist[i]);
		if (He3) {
			tmpMeas = new He3_Measurement(dataLoc + c_name + ".spc");
		} else {
			tmpMeas = new NaI_Measurement(dataLoc + c_name + ".Spc");
		}
		measCounts[i] = tmpMeas->GetTotalPulses();
		measTime[i] = tmpMeas->GetLiveTime();
		distData[i] = lexical_cast<double>(dist[i]);
		
		delete tmpMeas;
	}
	cpsData = (measCounts / measTime - bkg->GetCPS()) / activity;
	
	FitData();
}

DistResponse::DistResponse(const vector<boost::shared_ptr<BaseMeasurement>> distMeas, vector<double> dist, double activity, double activityUncertainty, const boost::shared_ptr<BaseMeasurement> bkg, bool curve_fit) : rand(chrono::system_clock::now().time_since_epoch().count()), measCounts(dist.size()), measTime(dist.size()), cpsData(dist.size()), distData(dist.size()), activity(activity), activityUncertainty(activityUncertainty), curve_fit(curve_fit) {
	if (distMeas.size() != dist.size()) {
		cout << "Warning: Number of distances and distance measurements are not equal!" << endl;
	}
	if (distMeas.size() == 0) {
		cout << "Error: No distance measurements!" << endl;
		return;
	}
	
	for (int i = 0; i < distMeas.size(); i++) {
		measCounts[i] = distMeas[i]->GetTotalPulses();
		measTime[i] = distMeas[i]->GetLiveTime();
		distData[i] = dist[i];
	}
	cpsData = (measCounts / measTime - bkg->GetCPS()) / activity;
	FitData();
}

DistResponse::DistResponse() : rand(chrono::system_clock::now().time_since_epoch().count()), p1(0.0), p2(0.0) {
	
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

double DistResponse::operator()(const double &dist) {
	return p1 / (4.0 * pi * pow(dist, p2));
}

Eigen::ArrayXd DistResponse::operator()(const Eigen::ArrayXd &dist) {
	return p1 / (4.0 * pi * dist.pow(p2));
}
