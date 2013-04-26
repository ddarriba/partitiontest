/*
 * Observable.cc
 *
 *  Created on: 30/05/2012
 *      Author: diego
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

void Observable::notify_observers(ObservableInfo * info, string * message) {
	for (unsigned int i = 0; i < observers.size(); i++) {
		((Observer *) observers[i])->update(*info);
	}
	delete info;
}

void Observable::notify_observers(MessageType type, unsigned int p_index,
		Model * model, time_t time, unsigned int current_index,
		unsigned int max_index, string * message) {
	notify_observers(
			new ObservableInfo(type, p_index, model, time, current_index,
					max_index, message));
}

void Observable::notify_observers(MessageType type, unsigned int p_index,
		time_t time, unsigned int current_index, unsigned int max_index, string * message) {
	ObservableInfo * oi = new ObservableInfo(type, p_index, (Model *) NULL,
			time, current_index, max_index, message);
	notify_observers(oi);
}

void Observable::notify_observers(MessageType type, unsigned int p_index,
		ModelSet * modelset, time_t time, unsigned int current_index,
		unsigned int max_index, string * message) {
	ObservableInfo * info = new ObservableInfo(type, p_index, modelset, time,
			current_index, max_index, message);
	notify_observers(info);
}

void Observable::chainObservers(Observable & item) {
	for (unsigned int i = 0; i < observers.size(); i++)
		item.attach(this->observers[i]);
}

} /* namespace partest */
