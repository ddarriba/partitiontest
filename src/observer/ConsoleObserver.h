/*
 * ConsoleObserver.hpp
 *
 *  Created on: 30/05/2012
 *      Author: diego
 */

#ifndef CONSOLEOBSERVER_H_
#define CONSOLEOBSERVER_H_

#include "Observer.h"
#include "../options/ParTestOptions.h"
#include <time.h>
//#include <boost/thread/mutex.hpp>

namespace partest {

class ConsoleObserver: public Observer {
  public:
    ConsoleObserver();
    virtual ~ConsoleObserver();
    virtual void update(const ObservableInfo & info,
        ParTestOptions * run_instance = NULL);
    virtual void update(string name, unsigned int current_index, unsigned int max_index);
  private:
    unsigned int numberOfTasks;
    unsigned int modelDigitsCount;
    time_t modelsetInitTime;
    time_t initTime;
    time_t endTime;
   // boost::mutex m_io_monitor;
    int currentScheme;
    unsigned int modelsetIndex;
    unsigned int modelsetCount;
    unsigned int modelIndex;
};

} /* namespace partest */
#endif /* CONSOLEOBSERVER_H_ */
