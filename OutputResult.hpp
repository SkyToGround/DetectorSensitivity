//
//  OutputResult.hpp
//

#ifndef OutputResult_hpp
#define OutputResult_hpp

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <eigen3/Eigen/Dense>
#include <sstream>
#include "ConsolePrint.hpp"

class OutputResult {
public:
	enum class OutputType {SCREEN, JSON_FILE};
	OutputResult(OutputType outType, std::string fileName = "");
	~OutputResult();
	void SetCoutKeyWidth(unsigned int width);
	
	void StartResult();
	void EndResult();
	template<typename T>
	void Write(const std::string &key, const T &value) {
		switch (outType) {
			case OutputType::SCREEN:
				CoutWrite(key, value);
				break;
			case OutputType::JSON_FILE:
				JsonWrite(key, value);
				break;
		}
	};
private:
	OutputType outType;
	unsigned int keyWidth;
	std::ofstream outStream;
	std::ostringstream tOut;
	bool firstJsonItm;
	bool firstDataItm;
	
	void CoutStartResult();
	void JsonStartResult();
	void CoutEndResult();
	void JsonEndResult();
	
	void CoutWrite(std::string key, std::string value);
	void JsonWrite(std::string key, std::string value);
	
	void CoutWrite(std::string key, double value);
	void JsonWrite(std::string key, double value);
	
	void CoutWrite(std::string key, int value);
	void JsonWrite(std::string key, int value);
	
	void CoutWrite(std::string key, bool value);
	void JsonWrite(std::string key, bool value);
	
	void CoutWrite(std::string key, std::vector<double> value);
	void JsonWrite(std::string key, std::vector<double> value);
	
	void CoutWrite(std::string key, Eigen::ArrayXd value);
	void JsonWrite(std::string key, Eigen::ArrayXd value);
	
	void CoutWrite(std::string key, std::vector<unsigned int> value);
	void JsonWrite(std::string key, std::vector<unsigned int> value);
	
	void CoutWrite(std::string key, std::vector<std::pair<double, double>> value);
	void JsonWrite(std::string key, std::vector<std::pair<double, double>> value);
	
	void JsonCheckFirst();
};

#endif /* OutputResult_hpp */
