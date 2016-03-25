//
//  main.cpp

#include <iostream>
#include "Detector.h"
#include <string>
#include <vector>
#include "Response.h"
#include "CalcCoordinator.h"
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <exception>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

using namespace std;

class

AngularResponse parseAngularData(std::vector<std::string> angularData, double bkg, bool &uncertainty);
DistResponse parseDistData(std::vector<std::string> distData, double bkg, double activity, double activity_uncertinaty, bool curve_fit, bool &uncertainty);
void parseActivity(std::string actString, double &activity, double &actUncert, bool &uncertainty);
void parseBkg(std::string bkgString, double &pulses, double &livetime, bool &uncertainty);
void parseActVsTime(std::string inString, double &start_time, double &stop_time, unsigned int &steps);

int main(int argc, const char * argv[])
{
	po::options_description generic("Generic arguments");
	generic.add_options()
	("help", "Show this help message")
	("config_file", po::value<std::string>(), "Retrive input parameters from a config file.");
	po::options_description conf("Configuration");
	conf.add_options()
	("fph", po::value<double>()->default_value(1.0), "Acceptable number of false positives per hour. Must be > 0. Default value is 1.")
	("beta", po::value<double>()->default_value(1/20), "Probability of false negative. Must be a value > 0 and < 1. Default value is 0.05.")
	("uncertainty", po::value<bool>()->default_value(false), "Perform uncertainty calculation. Possible values are 'true' and 'false'. Default is 'false'.")
	("uncertainty_loops", po::value<int>()->default_value(100), "Number of loops in uncertainty calculation. Default value is 100.")
	("calc_type", po::value<std::vector<std::string> >()->multitoken(), "Type of calculation to perform. Possible values are: 'mean', 'best', 'worst' and 'list'. Accepts multiple options. Defaults to 'mean'.")
	("output", po::value<string>(), "If argument is given, will output to given file in JSON-format. If not; outputs results to screen.")
	("distance", po::value<std::vector<double> >()->multitoken(), "The distance at which the source is placed in metres. Must be > 0. Accepts multiple values. Default value is 10.")
	("velocity", po::value<std::vector<double> >()->multitoken(), "The velocity of the detector in metres per second. Must be > 0. Accepts multiple values. Default value is 8.33 (30 kph).")
	("background", po::value<std::vector<double> >()->multitoken(), "The simulated background (in cps). Must be > 0. Accepts multiple values.")
	("cal_act", po::value<std::string>(), "The activity of the source used for calibration. Can be in the form of a single value or \"activity:uncertainty\". The uncertainty must be given in percent.")
	("cal_bkg", po::value<std::string>(), "The background of the detector in the environment where calibration measurements were done. Can be in the form of a single cps value or \"pulses:livetime\".")
	("dist_response", po::value<std::vector<std::string> >()->multitoken(), "Detector response input. Can be in one of two formats: \"distance:pulses:time\" or \"distance:cps\" where \"distance\" is the distance between source and detector in metres, \"pulses\" is the number of pulses reigstered by the detector, \"time\" is the live time of the measurement, and \"cps\" is the number of pulses per second. Accepts multiple values. See manual and examples for more information.")
	("dist_model", po::value<std::string>()->default_value("mean"), "How the distance response function is calculated. Possible values are \"least_sq\" and \"mean\". Defaults to \"mean\".")
	("ang_response", po::value<std::vector<std::string> >()->multitoken(), "Detector response input. Can be in one of two formats: \"angle:pulses:time\" or \"angle:cps\" where \"angle\" is in radians, \"pulses\" is the number of pulses reigstered by the detector, \"time\" is the live time of the measurement, and \"cps\" is the number of pulses per second. Accepts multiple values. Defaults to a angular response function of 1.0 over all angles; See manual and examples for more information.")
	("curve_limit", po::value<double>()->default_value(0.01), "When performing calculations of type 'best', 'mean' or 'worst'; for how long should the dynamic part of the intensity function be followed? See manual and examples for more information. Must be < 1.0 and > 0.0. Defaults to 0.01 (1% of peak value).")
	("mean_iters", po::value<unsigned int>()->default_value(100), "When performing calculations of type 'mean'; how many different time alignments should be tested in order to calculate the mean value? Defaults to 100.")
	("act_vs_time", po::value<std::string>(), "Create detection limit as a function of integration time data. Arguments must be of the form: \"start_time:stop_time:steps\". If this option is used, use of the \"output\" argument is strongly recommended. See manual and example files for more information.")
	("integration_time", po::value<std::vector<double> >()->multitoken(), "Set a fixed integration time and calculate minimum detectable activity from that parameter. Can take multiple values. Do not use together with the \"act_vs_time\" argument. Must be greater than 0.")
	("sim_iters", po::value<unsigned int>()->default_value(500000), "When performing calculations on list mode measurements, how many iterations should be used to calculate the probability of detecting a source. Default value is 20000.")
	;
	
	po::options_description cmd_line;
	cmd_line.add(generic).add(conf);
	
	po::options_description conf_file;
	conf_file.add(conf);
	
	po::variables_map vm;
	bool enableUncertainty = true;
	
	try {
	
		po::store(po::parse_command_line(argc, argv, cmd_line), vm);
		po::notify(vm);
	} catch (po::error e) {
		cout << "Unknown or invalid command line option. Try the \"--help\" option for a list of commands." << endl;
		
		return 0;
	}
	
	if (vm.count("config_file") > 0) {
		std::string fileName = vm["config_file"].as<std::string>();
		if (!boost::filesystem::exists(fileName)) {
			cout << "Unable to open config file, it does not exist." << endl;
			return 0;
		}
		ifstream ifs(fileName.c_str());
		if(!ifs) {
			cout << "Unable to open config file." << endl;
			return 0;
		} else {
			try {
				po::store(parse_config_file(ifs, conf_file), vm);
			} catch (boost::program_options::invalid_config_file_syntax e) {
				cout << "Unable to parse config file: \"" << e.what() << "\"" << endl;
				return 0;
			}
			po::notify(vm);
		}
	}
	
	if (vm.count("help") or argc == 1) {
		cout << cmd_line << "\n";
		return 0;
	}
	
	double fph = vm["fph"].as<double>();
	if (fph <= 0) {
		cout << "Number of false positives per hour (fph) must be more than 0." << endl;
		return 0;
	}
	
	double beta = vm["beta"].as<double>();
	if (beta <= 0 or beta >= 1.0) {
		cout << "Probability of false negative (beta) must fall within the range 0.0 < beta < 1.0." << endl;
		return 0;
	}
	
	std::vector<Detector::CalcType> cTypes;
	
	if (vm.count("calc_type") == 0) {
		cTypes.push_back(Detector::CalcType::MEAN);
	} else {
		std::vector<std::string> calcs = vm["calc_type"].as<std::vector<std::string>>();
		for (int y = 0; y < calcs.size(); y++) {
			if (calcs[y] == std::string("mean")) {
				cTypes.push_back(Detector::CalcType::MEAN);
			} else if (calcs[y] == std::string("worst")) {
				cTypes.push_back(Detector::CalcType::WORST);
			} else if (calcs[y] == std::string("best")) {
				cTypes.push_back(Detector::CalcType::BEST);
			} else if (calcs[y] == std::string("list")) {
				cTypes.push_back(Detector::CalcType::LIST_MODE);
			} else {
				cout << "Unknown calculation type \"" << calcs[y] << "\"." << endl;
				return 0;
			}
		}
		if (cTypes.size() == 0) {
			cout << "One or more type of calculation must be used." << endl;
			return 0;
		}
	}
	
	std::string output_file;
	if (vm.count("output") == 1) {
		output_file = vm["output"].as<std::string>();
	}
	
	std::vector<double> dists;
	if (vm.count("distance")) {
		dists = vm["distance"].as<vector<double> >();
	} else {
		dists.push_back(10.0);
	}
	for (int i = 0; i < dists.size(); i++) {
		if (dists[i] <= 0.0){
			cout << "Distances can not be negative." << endl;
			return 0;
		}
	}
	
	std::vector<double> velocities;
	if (vm.count("velocity")) {
		velocities = vm["velocity"].as<vector<double> >();
	} else {
		velocities.push_back(8.33);
	}
	for (int o = 0; o < velocities.size(); o++) {
		if (velocities[o] <= 0.0){
			cout << "Velocities can not be negative." << endl;
			return 0;
		}
	}
	
	double bkg_pulses, bkg_livetime;
	if (vm.count("cal_bkg") == 0) {
		cout << "Calibration background data is needed to run the simulation." << endl;
		return 0;
	} else {
		try {
			parseBkg(vm["cal_bkg"].as<std::string>(), bkg_pulses, bkg_livetime, enableUncertainty);
		} catch (std::runtime_error e) {
			cout << "Error when parsing calibration background data: \"" << e.what() << "\"" << endl;
			return 0;
		}
	}
	BkgResponse calBkg(bkg_pulses, bkg_livetime);
	
	double act, act_uncert;
	if (vm.count("cal_act") == 0) {
		cout << "Calibration source activity data is needed to run the simulation." << endl;
		return 0;
	} else {
		try {
			parseActivity(vm["cal_act"].as<std::string>(), act, act_uncert, enableUncertainty);
		} catch (std::runtime_error e) {
			cout << "Error when parsing activity data: \"" << e.what() << "\"" << endl;
			return 0;
		}
	}
	
	bool dist_fit_curve = false;
	if (vm.count("dist_model") != 0) {
		std::string distModel = vm["dist_model"].as<std::string>();
		if (std::string("mean") == distModel) {
			dist_fit_curve = false;
		} else if (std::string("least_sq") == distModel) {
			dist_fit_curve = true;
		} else {
			cout << "Error when reading distance model type. The model is unknown." << endl;
			return 0;
		}
	}
	
	DistResponse distResp;
	std::vector<std::string> distData;
	if (vm.count("dist_response") == 0) {
		cout << "Distance response data is needed to run the simulation." << endl;
		return 0;
	} else {
		distData	= vm["dist_response"].as<std::vector<std::string> >();
		try {
			distResp = parseDistData(distData, bkg_pulses / bkg_livetime, act, act_uncert, dist_fit_curve, enableUncertainty);
		} catch (std::runtime_error e) {
			cout << "Error when parsing distance data: \"" << e.what() << "\"" << endl;
			return 0;
		}
	}
	
	AngularResponse angResp;
	std::vector<std::string> angData;
	if (vm.count("ang_response") == 0) {
		
	} else {
		angData = vm["ang_response"].as<std::vector<std::string> >();
		try {
			angResp = parseAngularData(angData, bkg_pulses / bkg_livetime, enableUncertainty); //Fix me, bkg!
		} catch (std::runtime_error e) {
			cout << "Error when parsing angular data: \"" << e.what() << "\"" << endl;
			return 0;
		}
	}
	
	std::vector<double> background;
	if (vm.count("background") == 0) {
		cout << "One or more background count rates are needed to run the simulation." << endl;
		return 0;
	} else {
		background = vm["background"].as<std::vector<double> >();
		for (int i = 0; i < background.size(); i++) {
			if (background[i] <= 0.0) {
				cout << "Simulated background can not be equal to or less than zero." << endl;
				return 0;
			}
		}
	}
	
	bool uncert_calc = vm["uncertainty"].as<bool>();
	if (uncert_calc and !enableUncertainty) {
		cout << "Unable to calculate uncertainty as data needed to do so is not provided." << endl;
		return 0;
	}
	int uncertainty_loops = vm["uncertainty_loops"].as<int>();
	if (uncert_calc and uncertainty_loops < 2) {
		cout << "The number of uncertainty loops can not be less than 2." << endl;
		return 0;
	}
	if (!uncert_calc) {
		uncertainty_loops = 0;
	}
	
	double curve_limit = vm["curve_limit"].as<double>();
	if (curve_limit <= 0 or curve_limit >= 1.0) {
		cout << "The curve limit must be > 0 and < 1.0." << endl;
		return 0;
	}
	
	unsigned int mean_iters = vm["mean_iters"].as<unsigned int>();
	if (mean_iters == 0) {
		cout << "The number of iterations for calculating the mean must be greater than 0." << endl;
		return 0;
	}
	bool plotCalc = false;
	double start_time = 0, stop_time = 0;
	unsigned int steps = 0;
	if (vm.count("act_vs_time") == 1) {
		plotCalc = true;
		std::string actVsTimeData = vm["act_vs_time"].as<std::string>();
		try {
			parseActVsTime(actVsTimeData, start_time, stop_time, steps);
		} catch (std::runtime_error e) {
			cout << "Error when parsing activity vs time data: \"" << e.what() << "\"" << endl;
			return 0;
		}
	}
	
	bool fixedInt = false;
	std::vector<double> fixedIntTimes;
	if (vm.count("integration_time") > 0 and not plotCalc) {
		fixedInt = true;
		fixedIntTimes = vm["integration_time"].as<std::vector<double>>();
		for (int j = 0; j < fixedIntTimes.size(); j++) {
			if (fixedIntTimes[j] <= 0) {
				std::cout << "Integration time must be greater than 0." << std::endl;
				return 0;
			}
		}
	} else if (vm.count("integration_time") == 1 and plotCalc) {
		cout << "Can not use argument \"integration_time\" and \"act_vs_time\". Ignoring \"integration_time\"." << endl;
	}
	
	unsigned int sim_iters = vm["sim_iters"].as<unsigned int>();
	if (mean_iters == 0) {
		cout << "The number of simulation iterations must be greater than 0." << endl;
		return 0;
	}
	
	OutputResult::OutputType outDev = OutputResult::OutputType::SCREEN;
	if (output_file.size() != 0) {
		outDev = OutputResult::OutputType::JSON_FILE;
	}
	
	
	CalcCoordinator calc(outDev, output_file);
	Detector det(calBkg, distResp, angResp, curve_limit, mean_iters, sim_iters);
	
	for (int a = 0; a < dists.size(); a++) {
		det.SetDistance(dists[a]);
		for (int b = 0; b < velocities.size(); b++) {
			det.SetVelocity(velocities[b]);
			for (int c = 0; c < background.size(); c++) {
				det.SetSimBkg(background[c]);
				for (int d = 0; d < cTypes.size(); d++) {
					if (plotCalc) {
						calc.AddPlotCalculation(det, cTypes[d], fph, beta, start_time, stop_time, steps);
					} else if (fixedInt) {
						for (int e = 0; e < fixedIntTimes.size(); e++) {
							calc.AddFixedCalculation(det, cTypes[d], fph, beta, fixedIntTimes[e]);
						}
					} else {
						calc.AddCalculation(det, cTypes[d], fph, beta, uncertainty_loops);
					}
				}
			}
		}
	}
	
	calc.RunCalculations();
	return 0;
}

void parseActVsTime(std::string inString, double &start_time, double &stop_time, unsigned int &steps) {
	boost::regex expr{"^(\\d+\\.?\\d*):(\\d+\\.?\\d*):(\\d+)$"};
	boost::smatch res;
	bool match = boost::regex_match(inString, res, expr);
	if (!match) {
		throw std::runtime_error(std::string("Incorrectly formated activity vs time parameters."));
	}
	start_time = lexical_cast<double>(res[1]);
	stop_time = lexical_cast<double>(res[2]);
	steps = lexical_cast<unsigned int>(res[3]);
	if (start_time >= stop_time) {
		throw std::runtime_error(std::string("Start time can not be equal to or greater than stop time."));
	}
}

void parseActivity(std::string actString, double &activity, double &actUncert, bool &uncertainty) {
	boost::regex expr{"^(\\d+\\.?\\d*)(?::(\\d+\\.?\\d*))?$"};
	boost::smatch res;
	bool match = boost::regex_match(actString, res, expr);
	if (!match) {
		throw std::runtime_error(std::string("Incorrectly formated calibration source activity data."));
	}
	activity = lexical_cast<double>(res[1]);
	if (activity <= 0) {
		throw std::runtime_error(std::string("Calibration source activity can not be less than or equal to zero."));
	}
	if (res[2].length() != 0) {
		actUncert = lexical_cast<double>(res[2]);
		if (actUncert < 0) {
			throw std::runtime_error(std::string("Calibration source activity uncertainty can not be less than zero."));
		}
	} else {
		actUncert = 0.0;
		uncertainty = false;
	}
}

void parseBkg(std::string bkgString, double &pulses, double &livetime, bool &uncertainty) {
	boost::regex expr{"^(\\d+\\.?\\d*)(?::(\\d+\\.?\\d*))?$"};
	boost::smatch res;
	bool match = boost::regex_match(bkgString, res, expr);
	if (!match) {
		throw std::runtime_error(std::string("Incorrectly formated calibration background data."));
	}
	pulses = lexical_cast<double>(res[1]);
	if (pulses < 0) {
		throw std::runtime_error(std::string("Background can not be less than zero."));
	}
	if (res[2].length() != 0) {
		livetime = lexical_cast<double>(res[2]);
		if (livetime <= 0) {
			throw std::runtime_error(std::string("Background live time can not be less than or equal to zero."));
		}
	} else {
		livetime = 1.0;
		uncertainty = false;
	}
}


AngularResponse parseAngularData(std::vector<std::string> angularData, double bkg, bool &uncertainty) {
	boost::regex expr{"^(\\d+\\.?\\d*):(\\d+\\.?\\d*)(?::(\\d+\\.?\\d*))?$"};
	boost::smatch res;
	std::vector<double> livetime;
	std::vector<double> pulses;
	std::vector<double> angle;
	for (int y = 0; y < angularData.size(); y++) {
		bool match = boost::regex_match(angularData[y], res, expr);
		if (!match) {
			throw std::runtime_error(std::string("Incorrectly formated angular data."));
		}
		angle.push_back(lexical_cast<double>(res[1]));
		pulses.push_back(lexical_cast<double>(res[2]));
		if (res[3].length() != 0) {
			livetime.push_back(lexical_cast<double>(res[3]));
		} else {
			livetime.push_back(1.0);
			uncertainty = false;
		}
		if (angle[angle.size() - 1] < 0 or angle[angle.size() - 1] > (pi / 1.99)) { //Fix me: should be 2 but using 1.99 for UI reasons
			throw std::runtime_error(std::string("Angular value does not fall within the range 0.0 to pi / 2."));
		}
		if (pulses[pulses.size() - 1] <= 0) {
			throw std::runtime_error(std::string("Number of pulses can not be zero or negative."));
		}
		
		if (livetime[livetime.size() - 1] <= 0) {
			throw std::runtime_error(std::string("Live time can not be zero or negative."));
		}
	}
	return AngularResponse(pulses, livetime, angle, bkg);
}

DistResponse parseDistData(std::vector<std::string> distData, double bkg, double activity, double activity_uncertinaty, bool curve_fit, bool &uncertainty) {
	boost::regex expr{"^(\\d+\\.?\\d*):(\\d+\\.?\\d*)(?::(\\d+\\.?\\d*))?$"};
	boost::smatch res;
	std::vector<double> livetime;
	std::vector<double> pulses;
	std::vector<double> distance;
	for (int y = 0; y < distData.size(); y++) {
		bool match = boost::regex_match(distData[y], res, expr);
		if (!match) {
			throw std::runtime_error(std::string("Incorrectly formated distance data."));
		}
		distance.push_back(lexical_cast<double>(res[1]));
		pulses.push_back(lexical_cast<double>(res[2]));
		if (res[3].length() != 0) {
			livetime.push_back(lexical_cast<double>(res[3]));
		} else {
			livetime.push_back(1.0);
			uncertainty = false;
		}
		if (distance[distance.size() - 1] <= 0) {
			throw std::runtime_error(std::string("Distance value can not be less than or equal to zero."));
		}
		if (pulses[pulses.size() - 1] <= 0) {
			throw std::runtime_error(std::string("Number of pulses can not be zero or negative."));
		}
		
		if (livetime[livetime.size() - 1] <= 0) {
			throw std::runtime_error(std::string("Live time can not be zero or negative."));
		}
	}
	return DistResponse(pulses, livetime, distance, bkg, activity, activity_uncertinaty, curve_fit);
}