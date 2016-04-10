//
//  concurrent_queue.h

#ifndef _CONCURRENT_QUEUE_H_
#define _CONCURRENT_QUEUE_H_

#include <queue>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/xtime.hpp>

template<typename Data>
class concurrent_queue
{
private:
	std::queue<Data> the_queue;
	mutable boost::mutex the_mutex;
	boost::condition_variable the_condition_variable;
public:
	void push(Data const& data)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		the_queue.push(data);
		lock.unlock();
		the_condition_variable.notify_one();
	}
	
	unsigned int size()
	{
		unsigned int retVal;
		boost::mutex::scoped_lock(the_mutex);
		retVal = the_queue.size();
		return retVal;
	}
	
	bool empty() const
	{
		boost::mutex::scoped_lock lock(the_mutex);
		return the_queue.empty();
	}
	
	bool try_pop(Data& popped_value)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		if(the_queue.empty())
		{
			return false;
		}
		
		popped_value=the_queue.front();
		the_queue.pop();
		return true;
	}
	
	void wait_and_pop(Data& popped_value)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		while(the_queue.empty())
		{
			the_condition_variable.wait(lock);
		}
		
		popped_value=the_queue.front();
		the_queue.pop();
	}
	
	bool time_out_pop(Data& popped_value)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		boost::xtime td;
		xtime_get(&td, boost::TIME_UTC_);
		td.nsec += 200000000;
		if (the_queue.empty())
		{
			the_condition_variable.timed_wait(lock, td);
		}
		
		if (the_queue.empty())
		{
			return false;
		}
		
		popped_value=the_queue.front();
		the_queue.pop();
		return true;
	}
	
};

#endif