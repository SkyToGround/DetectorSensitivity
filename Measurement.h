//
//  BaseMeasurement.h

#ifndef __NeutronDetectorSim__BaseMeasurement__
#define __NeutronDetectorSim__BaseMeasurement__

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <random>
#include <chrono>


using namespace Eigen;

using namespace std;

class BaseMeasurement {
public:
	BaseMeasurement();
	virtual ~BaseMeasurement();
	virtual double GetUncertainty() = 0;
	virtual double GetCPS() = 0;
	virtual double GetTotalPulses() = 0;
	virtual double GetLiveTime() = 0;
	virtual double GetRandomizedCPS() = 0;
};

class NaI_Measurement : public BaseMeasurement {
public:
	NaI_Measurement(string fileName);
	NaI_Measurement(const NaI_Measurement &setMeas);
	double GetCPS();
	double GetUncertainty();
	double GetTotalPulses() {return totalPulses;};
	double GetLiveTime() {return liveTime;};
	double GetRandomizedCPS();
private:
#pragma pack(2)
	struct headerData {
		short INFTYPE;//1
		short file_type;//2
		short reserved[2];//3,4
		short acq_ptr;//5
		short sample_desc_ptr;//6
		short det_desc_ptr;//7
		short EBR_desc_ptr;//8
		short an1_par_ptr;//9
		short an2_par_ptr;//10
		short an3_par_ptr;//11
		short an4_par_ptr;//12
		short abs_corr_desc_ptr;//13
		short IEQ_desc_ptr;//14
		short geo_desc_ptr;//15
		short mpc_desc_ptr;//16
		short cal_desc_ptr;//17
		short cal1_data_ptr;//18
		short cal2_data_ptr;//19
		short first[11];
		
		short specPtr;
		short specCtr;
		short channels;
		short something;
		float second;
		double third;
		short fourth[5];
		float realtime;
		float livetime;
	};
	
	struct calData {
		short above_knee_eff_type;//1
		short below_knee_eff_type;//2
		short eff_pairs;//3
		short channels;//4
		float knee;//5
		float sig_abv_knee;//7
		float sig_bel_knee;//9
		float en_A;//11
		float en_B;//13
		float en_C;//15
	};
#pragma pack()
	double totalPulses;
	double liveTime;
	double cps;
};

class He3_Measurement : public BaseMeasurement {
public:
	He3_Measurement(string fileName);
	double GetCPS() {return cps;};
	double GetUncertainty() {return 1.0;};
	double GetTotalPulses() {return totalPulses;};
	double GetLiveTime() {return liveTime;};
	double GetRandomizedCPS();
private:
#pragma pack(1)
	typedef struct {
		unsigned long	record_number;
		unsigned long	line_number;
		double			UTC_time;
		double			X;
		double			Y;
		double			Z;
		double			northing;
		double			easting;
		float			PDOP;
		long			DGPS;
		float			ralt;
		float			live_time;
		float			roi[10];
		short int		gain;
		short int		peak;
		unsigned short	spec[512];
	} SPCdata;
#pragma pack(0)
	void ExtractData(char *data, long long int length, ArrayXd &pulses, ArrayXd &samples);
	double totalPulses;
	double liveTime;
	double totalSamples;
	double cps;
};

#endif /* defined(__NeutronDetectorSim__BaseMeasurement__) */
