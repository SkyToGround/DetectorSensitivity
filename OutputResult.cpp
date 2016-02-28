//
//  OutputResult.cpp
//  MobileDetectorSim
//
//  Created by Jonas Nilsson on 2016-02-28.
//  Copyright © 2016 Jonas Nilsson. All rights reserved.
//

#include "OutputResult.hpp"

OutputResult::OutputResult(OutputType outType, std::string fileName) : outType(outType), outStream(), firstJsonItm(false), keyWidth(15) {
	if (outType == OutputType::JSON_FILE) {
		outStream.open(fileName);
		//Fix me: check if stream is good
		outStream << "[\n";
	}
}

OutputResult::~OutputResult() {
	if (outType == OutputType::JSON_FILE) {
		outStream << "\n]\n";
		outStream.close();
	}
}

void OutputResult::SetCoutKeyWidth(unsigned int width) {
	keyWidth = width;
}

void OutputResult::StartResult() {
	switch (outType) {
		case OutputType::SCREEN:
			CoutStartResult();
			break;
		case OutputType::JSON_FILE:
			JsonStartResult();
			break;
	}
}

void OutputResult::CoutStartResult() {
	std::cout << "\n------------------";
}

void OutputResult::JsonStartResult() {
	if (not firstJsonItm) {
		firstJsonItm = true;
	} else {
		outStream << ",";
	}
	outStream << "\n\t{";
	firstDataItm = false;
}

void OutputResult::EndResult() {
	switch (outType) {
		case OutputType::SCREEN:
			CoutEndResult();
			break;
		case OutputType::JSON_FILE:
			JsonEndResult();
			break;
	}
}

void OutputResult::CoutEndResult() {
	std::cout << "\n------------------\n";
	std::cout.flush();
}

void OutputResult::JsonEndResult() {
	outStream << "\n\t}";
	outStream.flush();
}

void OutputResult::JsonCheckFirst() {
	if (not firstDataItm) {
		firstDataItm = true;
	} else {
		outStream << ",";
	}
}

void OutputResult::CoutWrite(std::string key, std::string value) {
	std::cout << "\n" << std::setw(keyWidth) << key << ": " << value;
}

void OutputResult::JsonWrite(std::string key, std::string value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": \"" << value << "\"";
}

void OutputResult::CoutWrite(std::string key, double value) {
	std::cout << "\n" << std::setw(keyWidth) << key << ": " << value;
}

void OutputResult::JsonWrite(std::string key, double value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": " << value;
}

void OutputResult::CoutWrite(std::string key, int value) {
	std::cout << "\n" << std::setw(keyWidth) << key << ": " << value;
}

void OutputResult::JsonWrite(std::string key, int value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": " << value;
}

void OutputResult::CoutWrite(std::string key, bool value) {
	if (value) {
		std::cout << "\n" << std::setw(keyWidth) << key << ": true";
	} else {
		std::cout << "\n" << std::setw(keyWidth) << key << ": false";
	}
	
}

void OutputResult::JsonWrite(std::string key, bool value) {
	JsonCheckFirst();
	if (value) {
		outStream << "\n\t\t\"" << key << "\": true";
	} else {
		outStream << "\n\t\t\"" << key << "\": false";
	}
}

void OutputResult::CoutWrite(std::string key, std::vector<double> value) {
	std::cout << "\n<<<<  " << key << "  >>>>\n";
	for (int i = 0; i < value.size(); i++) {
		std::cout << i + 1 << " : " << value[i] << "\n";
	}
	std::cout << "--  --  --  --";
}

void OutputResult::JsonWrite(std::string key, std::vector<double> value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": [" << value[0];
	for (int i = 1; i < value.size(); i++) {
		outStream << ", " << value[i];
	}
	outStream << "]";
}

void OutputResult::CoutWrite(std::string key, Eigen::ArrayXd value) {
	std::cout << "\n<<<<  " << key << "  >>>>\n";
	for (int i = 0; i < value.size(); i++) {
		std::cout << i + 1 << " : " << value[i] << "\n";
	}
	std::cout << "--  --  --  --";
}

void OutputResult::JsonWrite(std::string key, Eigen::ArrayXd value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": [" << value[0];
	for (int i = 1; i < value.size(); i++) {
		outStream << ", " << value[i];
	}
	outStream << "]";
}

void OutputResult::CoutWrite(std::string key, std::vector<unsigned int> value) {
	std::cout << "\n<<<<  " << key << "  >>>>\n";
	for (int i = 0; i < value.size(); i++) {
		std::cout << i + 1 << " : " << value[i] << "\n";
	}
	std::cout << "--  --  --  --";
}

void OutputResult::JsonWrite(std::string key, std::vector<unsigned int> value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": [" << value[0];
	for (int i = 1; i < value.size(); i++) {
		outStream << ", " << value[i];
	}
	outStream << "]";
}

void OutputResult::CoutWrite(std::string key, std::vector<std::pair<double, double>> value) {
	std::cout << "\n<<<<  " << key << "  >>>>\n";
	for (int i = 0; i < value.size(); i++) {
		std::cout << value[i].first << " : " << value[i].second << "\n";
	}
	std::cout << "--  --  --  --";
}

void OutputResult::JsonWrite(std::string key, std::vector<std::pair<double, double>> value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": [[" << value[0].first << ", " << value[0].second << "]";
	for (int i = 1; i < value.size(); i++) {
		outStream << ", [" << value[i].first << ", " << value[i].second << "]";
	}
	outStream << "]";
}

//void OutputResult::WriteB(const std::string &key, const bool &value) {
//	switch (outType) {
//		case OutputType::SCREEN:
//			CoutWriteB(key, value);
//			break;
//		case OutputType::JSON_FILE:
//			JsonWriteB(key, value);
//			break;
//	}
//};
