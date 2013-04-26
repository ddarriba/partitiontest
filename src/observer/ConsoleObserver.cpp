/*
 * ConsoleObserver.cc
 *
 *  Created on: 30/05/2012
 *      Author: diego
 */

#include "ConsoleObserver.h"
#include "../util/Utilities.h"
#include "../util/GlobalDefs.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>

namespace partest {

ConsoleObserver::ConsoleObserver() {
	init_time = time(NULL);
	current_partition = -1;
	modelsetCount = modelsetIndex = model_digits_count = modelIndex =
			number_of_tasks = 0;
}

ConsoleObserver::~ConsoleObserver() {
}

void ConsoleObserver::update(string name, unsigned int current_index,
		unsigned int max_index) {
	cout << setw(11) << Utilities::timeToString(time(NULL) - init_time)
			<< " - - - Modelset " << current_index + 1 << "/" << max_index
			<< "  " << name << endl;
}
void ConsoleObserver::update(const ObservableInfo & info,
		ParTestOptions * run_instance) {
	//boost::mutex::scoped_lock lock(m_io_monitor);

	switch (info.type) {
	case MT_SINGLE_INIT:
#ifdef DEBUG
		cout << "[TRACE] [" << modelIndex+1 << "/" << info.max_index << "]\tOptimizing " << info.model->getName() << "..." << endl;
#endif
		break;
	case MT_SINGLE_END:
		if (modelIndex == 0)
			cout << "INIT MODEL" << endl;
		cout << "[" << ++modelIndex << "/" << info.max_index << "]\tOptimized "
				<< setw(15) << info.model->getName() << info.model->getLnL()
				<< endl;
		break;
	case MT_MODELSET_INIT:
		cout << "INIT MODELSET " << info.current_index << endl;
		modelIndex = 0;
		break;
	case MT_MODELSET_END:
		cout << "END MODELSET" << endl;
		break;
	case MT_PARTITION_INIT:
		cout << "INIT PARTITION " << *info.message << endl;
		modelIndex = 0;
		break;
	case MT_PARTITION_END:
		cout << "END PARTITION " << *info.message << endl;
		break;
	}
}

} /* namespace partest */
