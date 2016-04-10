//
//  OutputResult.cpp
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
	tOut.str(std::string(""));
	tOut << "\n------------------";
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
	tOut << "\n------------------\n";
	PRINT(tOut.str());
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
	tOut << "\n" << std::setw(keyWidth) << key << ": " << value;
}

void OutputResult::JsonWrite(std::string key, std::string value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": \"" << value << "\"";
}

void OutputResult::CoutWrite(std::string key, double value) {
	tOut << "\n" << std::setw(keyWidth) << key << ": " << value;
}

void OutputResult::JsonWrite(std::string key, double value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": " << value;
}

void OutputResult::CoutWrite(std::string key, int value) {
	tOut << "\n" << std::setw(keyWidth) << key << ": " << value;
}

void OutputResult::JsonWrite(std::string key, int value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": " << value;
}

void OutputResult::CoutWrite(std::string key, bool value) {
	if (value) {
		tOut << "\n" << std::setw(keyWidth) << key << ": true";
	} else {
		tOut << "\n" << std::setw(keyWidth) << key << ": false";
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
	tOut << "\n<<<<  " << key << "  >>>>\n";
	for (int i = 0; i < value.size(); i++) {
		tOut << i + 1 << " : " << value[i] << "\n";
	}
	tOut << "--  --  --  --";
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
	tOut << "\n<<<<  " << key << "  >>>>\n";
	for (int i = 0; i < value.size(); i++) {
		tOut << i + 1 << " : " << value[i] << "\n";
	}
	tOut << "--  --  --  --";
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
	tOut << "\n<<<<  " << key << "  >>>>\n";
	for (int i = 0; i < value.size(); i++) {
		tOut << i + 1 << " : " << value[i] << "\n";
	}
	tOut << "--  --  --  --";
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
	tOut << "\n<<<<  " << key << "  >>>>\n";
	for (int i = 0; i < value.size(); i++) {
		tOut << value[i].first << " : " << value[i].second << "\n";
	}
	tOut << "--  --  --  --";
}

void OutputResult::JsonWrite(std::string key, std::vector<std::pair<double, double>> value) {
	JsonCheckFirst();
	outStream << "\n\t\t\"" << key << "\": [[" << value[0].first << ", " << value[0].second << "]";
	for (int i = 1; i < value.size(); i++) {
		outStream << ", [" << value[i].first << ", " << value[i].second << "]";
	}
	outStream << "]";
}
