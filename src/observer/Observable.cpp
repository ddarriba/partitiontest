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
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

#include "Observable.h"

namespace partest {

Observable::Observable() {
}

Observable::~Observable() {
}

void Observable::attach(Observer *myObserver) {
	observers.push_back(myObserver);
}

void Observable::detach(Observer *myObserver) {
	for (unsigned int i = 0; i < observers.size(); i++) {
		if (observers[i] == myObserver) {
			observers.erase(observers.begin() + i);
			return;
		}
	}
}

void Observable::notify_observers(ObservableInfo * info, string message) {
	for (unsigned int i = 0; i < observers.size(); i++) {
		((Observer *) observers[i])->update(*info);
	}
	delete info;
}

void Observable::notify_observers(MessageType type,
		time_t time, string message) {
	t_partitionElementId nullId;
	notify_observers(
			new ObservableInfo(type, nullId, (Model *) NULL, time, 0, 0, message));
}

void Observable::notify_observers(MessageType type,
		t_partitionElementId p_index, Model * model, time_t time,
		unsigned int current_index, unsigned int max_index, string message) {
	notify_observers(
			new ObservableInfo(type, p_index, model, time, current_index,
					max_index, message));
}

void Observable::notify_observers(MessageType type,
		t_partitionElementId p_index, time_t time, unsigned int current_index,
		unsigned int max_index, string message) {
	ObservableInfo * oi = new ObservableInfo(type, p_index, (Model *) NULL,
			time, current_index, max_index, message);
	notify_observers(oi);
}

void Observable::notify_observers(MessageType type,
		t_partitionElementId p_index, ModelSet * modelset, time_t time,
		unsigned int current_index, unsigned int max_index, string message) {
	ObservableInfo * info = new ObservableInfo(type, p_index, modelset, time,
			current_index, max_index, message);
	notify_observers(info);
}

void Observable::chainObservers(Observable & item) {
	for (unsigned int i = 0; i < observers.size(); i++)
		item.attach(this->observers[i]);
}

} /* namespace partest */
