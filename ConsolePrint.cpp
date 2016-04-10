//
//  ConsolePrint.cpp

#include "ConsolePrint.hpp"

ConsolePrint *ConsolePrint::inst_pointer = NULL;

ConsolePrint* ConsolePrint::Inst()
{
	if (NULL == inst_pointer) {
		inst_pointer = new ConsolePrint(false);
	}
	return inst_pointer;
}

void ConsolePrint::Close()
{
	if (NULL != inst_pointer) {
		inst_pointer->Closer();
	}
}

void ConsolePrint::Closer()
{
	print_thread->interrupt();
	print_thread->join();
	delete print_thread;
	delete this;
}

void ConsolePrint::Init(bool printToConsole)
{
	if (NULL == inst_pointer) {
		inst_pointer = new ConsolePrint(printToConsole);
	}
}

ConsolePrint::~ConsolePrint()
{
}

ConsolePrint::ConsolePrint(bool printToConsole) : printToConsole(printToConsole) {
	print_thread = new boost::thread(&ConsolePrint::thread_function, this);
}

void ConsolePrint::Print(std::string message) {
	print_strct c_struct;
	
	c_struct.overridePrint = true;
	c_struct.message = message;
	print_queue.push(c_struct);
}

void ConsolePrint::ThreadPrint(std::string message, bool overridePrint) {
	print_strct c_struct;
	
	c_struct.overridePrint = overridePrint;
	c_struct.message = message + std::string("(Thread: ") + boost::lexical_cast<std::string>(boost::this_thread::get_id()) + std::string(")");
	print_queue.push(c_struct);
}

void ConsolePrint::thread_function()
{
	if (printToConsole) {
		std::cout << "Entering verbose mode." << std::endl;
	}
	print_strct c_msg;
	try {
		while (true) {
			boost::this_thread::interruption_point();
			if (print_queue.time_out_pop(c_msg)) {
				if (printToConsole or c_msg.overridePrint) {
					std::cout << c_msg.message << std::endl;
				}
				
				while (!print_queue.empty()) {
					print_queue.wait_and_pop(c_msg);
					if (printToConsole or c_msg.overridePrint) {
						std::cout << c_msg.message << std::endl;
					}
				}
			}
		}
	} catch (boost::thread_interrupted&) {
		if (printToConsole) {
			std::cout << "ConsolePrint::thread_function((): Recieved thread interrupt." << std::endl;
		}
	}
	if (printToConsole) {
		std::cout << "Ending verbose mode." << std::endl;
	}
}
