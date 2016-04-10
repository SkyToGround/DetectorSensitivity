//
//  ConsolePrint.hpp

#ifndef ConsolePrint_hpp
#define ConsolePrint_hpp

#include <string>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include "concurrent_queue.h"

#define PRINT ConsolePrint::Inst()->Print
#define PRINT_T ConsolePrint::Inst()->ThreadPrint

class ConsolePrint {
public:
	static void Init(bool printToConsole);
	static void Close();
	static ConsolePrint* Inst();
	void ThreadPrint(std::string message, bool overridePrint = false);
	void Print(std::string message);
	~ConsolePrint();
protected:
	ConsolePrint(bool printToConsole);
	ConsolePrint(ConsolePrint const&){};
private:
	struct print_strct {
		std::string message;
		bool overridePrint;
	};
	
	bool printToConsole;
	
	void Closer();
	void thread_function();
	boost::thread *print_thread;
	static ConsolePrint *inst_pointer;
	concurrent_queue<print_strct> print_queue;
};

//void PRINT(std::string message) {
//	ConsolePrint::Inst()->Print(message);
//}
//
//void PRINT_T(std::string message, bool overridePrint = false) {
//	ConsolePrint::Inst()->ThreadPrint(message, overridePrint);
//}

#endif /* ConsolePrint_hpp */
