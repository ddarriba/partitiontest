/*  PartitionTest, fast selection of the best fit partitioning scheme for
 *  multi-gene data sets.
 *  Copyright May 2013 by Diego Darriba
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other enquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file ConsoleObserver.h
 */

#ifndef CONSOLEOBSERVER_H_
#define CONSOLEOBSERVER_H_

#include "Observer.h"
#include "options/ParTestOptions.h"
#include <time.h>

namespace partest {

/**
 * @brief Observer that writes the output in the command console.
 */
class ConsoleObserver: public Observer {
public:
	ConsoleObserver();
	virtual ~ConsoleObserver();
	virtual void update(const ObservableInfo & info,
			ParTestOptions * run_instance = NULL);
	virtual void update(string name, unsigned int current_index,
			unsigned int max_index);
private:
	unsigned int numberOfTasks; /** Internal counter of tasks within a group. */
	unsigned int modelDigitsCount; /** Number of bits of a model */
	time_t modelsetInitTime; /** Internal timestamp of the beginning of a model set computation */
	time_t initTime; /** Internal timestamp */
	time_t endTime; /** Internal timestamp */
	int currentScheme; /** Internal mark of the current partitioning scheme. */
	unsigned int modelsetIndex; /** Internal index within a model set. */
	unsigned int modelsetCount; /** Internal size of a model set. */
	unsigned int modelIndex; /** Internal index of a model. */
};

} /* namespace partest */
#endif /* CONSOLEOBSERVER_H_ */
