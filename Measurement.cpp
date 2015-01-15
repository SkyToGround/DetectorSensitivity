//
//  BaseMeasurement.cpp

#include "Measurement.h"

BaseMeasurement::BaseMeasurement() {
	
}

BaseMeasurement::~BaseMeasurement() {
	
}

NaI_Measurement::NaI_Measurement(string fileName) {
	fstream inStream(fileName.c_str(), ios::in);
	if (!inStream.is_open()) {
		cout << "File \"" << fileName << "\" not open!" << endl;
		exit(-1);
	}
	headerData header;
	inStream.read((char*)&header, sizeof(header));
	liveTime = header.livetime;
	int channels = header.channels;
	unsigned int *spectra = new unsigned int[channels];
	
	calData enCal;
	inStream.seekg((header.cal1_data_ptr - 1) * 128, ios::beg);
	inStream.read((char*)&enCal, sizeof(enCal));
	
	double low_en = 2800;
	double high_en = 8000;
	int low_ch = int((low_en - enCal.en_A) / enCal.en_B + 0.5);
	int high_ch = int((high_en - enCal.en_A) / enCal.en_B + 0.5);
	
	inStream.seekg((header.specPtr - 1) * 128, ios::beg);
	inStream.read((char*)spectra, sizeof(unsigned int) * channels);
	totalPulses = 0;
//	for (int i = 270; i < 769; i++) {
//		totalPulses = totalPulses + spectra[i]; // Fix me
//	}
	
	for (int i = low_ch; i < high_ch; i++) {
		totalPulses = totalPulses + spectra[i];
	}
	cps = totalPulses / liveTime;
	
	delete [] spectra;
	inStream.close();
}

double NaI_Measurement::GetCPS() {
	return cps;
}

double NaI_Measurement::GetUncertainty() {
	return sqrt(totalPulses) / liveTime;
}

double NaI_Measurement::GetRandomizedCPS() {
	mt19937 rand(chrono::system_clock::now().time_since_epoch().count());
	
	normal_distribution<double> pulseDist(totalPulses, sqrt(totalPulses));
	
	return pulseDist(rand) / liveTime;
}

double He3_Measurement::GetRandomizedCPS() {
	mt19937 rand(chrono::system_clock::now().time_since_epoch().count());
	
	normal_distribution<double> pulseDist(totalPulses, sqrt(totalPulses));
	
	return pulseDist(rand) / liveTime;
}

He3_Measurement::He3_Measurement(string fileName) {
	ifstream inStream(fileName.c_str(), ios::in | ios::binary);
	if (!inStream.is_open()) {
		cout << "File \"" << fileName << "\" not open!" << endl;
		exit(-1);
	}
	inStream.seekg(0, ios::end);
	int dataLen = inStream.tellg();
	char *inData = new char[dataLen];
	inStream.seekg(0, ios::beg);
	inStream.read(inData, dataLen);
	inStream.close();
	
	ArrayXd pulses, samples;
	ExtractData(inData, dataLen, pulses, samples);
	delete [] inData;
	
	samples = samples / 10.0;
	
	totalPulses = pulses.sum();
	totalSamples = samples.sum();
	liveTime = totalSamples;
	cps = totalPulses / totalSamples;
}

void He3_Measurement::ExtractData(char *data, long long int length, ArrayXd &pulses, ArrayXd &samples) {
	char hHead[16];
	int first, second;
	memcpy(&hHead, data, 16);
	if (!(sscanf(hHead, "%d %d", &first, &second) == 2)) {
		cout << "Unable to identify as *.spc-file." << endl;
		return;
	}
	unsigned int currentFilePosition = 0;
	if (sizeof(SPCdata) == second) {
		currentFilePosition = first;
	} else {
		cout << "Incorrect package length! It should be " << second << " but it is " << sizeof(SPCdata) << endl;
		return;
	}
	
	SPCdata tempData;
	vector<double> tempSamples;
	vector<double> tempPulses;
	for (unsigned int t = 0; currentFilePosition < length; t++, currentFilePosition = currentFilePosition + sizeof(SPCdata)) {
		memcpy(&tempData, data + currentFilePosition, sizeof(SPCdata));
		tempSamples.push_back(tempData.spec[509]);
		tempPulses.push_back(tempData.spec[510]);
	}
	pulses = ArrayXd(tempPulses.size());
	samples = ArrayXd(tempPulses.size());
	for (int i = 0; i < tempPulses.size(); i++) {
		pulses[i] = tempPulses[i];
		samples[i] = tempSamples[i];
	}
	return;
}