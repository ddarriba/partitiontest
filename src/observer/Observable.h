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

/**
 * @file Observable.h
 */

#ifndef OBSERVABLE_H_
#define OBSERVABLE_H_

#include <vector>
#include <string>
#include "Observer.h"
#include "model/Model.h"

namespace partest {

/**
 * @brief Implementation for observable classes.
 */
class Observable {
public:
	Observable();
	virtual ~Observable();

	/**
	 * @brief Attaches a new Observer to the instance.
	 *
	 * @param[in] myObserver The new Observer.
	 */
	void attach(Observer * myObserver);

	/**
	 * @brief Detaches a new Observer from the observers of this instance.
	 *
	 * @param[in] myObserver The Observer to be detached.
	 */
	void detach(Observer * myObserver);

	/**
	 * @brief Notifies something to all the Observers
	 *
	 * @param[in] info The message to be sent.
	 * @param[in] message Human readable message.
	 */
	void notify_observers(ObservableInfo * info, string message = "");

	void notify_observers(MessageType type,	time_t time, string message = "");

	/**
	 * @brief Notifies something to all the Observers
	 *
	 * @param[in] type Message type.
	 * @param[in] p_index Index with some meaning for the observers.
	 * @param[in] model Reference to a model (if appropriate).
	 * @param[in] time Timestamp with some meaning for the observers.
	 * @param[in] current_index Current index within the group.
	 * @param[in] max_index Maximum index within the group.
	 * @param[in] message Human readable message.
	 */
	void notify_observers(MessageType type, t_partitionElementId p_index, Model * model,
			time_t time , unsigned int current_index,
			unsigned int max_index, string message = "");

	/**
	 * @brief Notifies something to all the Observers
	 *
	 * @param[in] type Message type.
	 * @param[in] p_index Index with some meaning for the observers.
	 * @param[in] time Timestamp with some meaning for the observers.
	 * @param[in] current_index Current index within the group.
	 * @param[in] max_index Maximum index within the group.
	 * @param[in] message Human readable message.
	 */
	void notify_observers(MessageType type, t_partitionElementId p_index, time_t time,
			unsigned int current_index, unsigned int max_index,
			string message = "");

	/**
	 * @brief Notifies something to all the Observers
	 *
	 * @param[in] type Message type.
	 * @param[in] p_index Index with some meaning for the observers.
	 * @param[in] modelset Reference to a modelset (if appropriate).
	 * @param[in] time Timestamp with some meaning for the observers.
	 * @param[in] current_index Current index within the group.
	 * @param[in] max_index Maximum index within the group.
	 * @param[in] message Human readable message.
	 */
	void notify_observers(MessageType type, t_partitionElementId p_index,
			ModelSet * modelset, time_t time,
			unsigned int current_index, unsigned int max_index,
			string message = "");

protected:

	/**
	 * @brief Attaches the observers to another Observable instance.
	 *
	 * @param item The Observable to attach the observers from this instance.
	 */
	void chainObservers(Observable & item);

private:
	vector<Observer *> observers; /** Set of Observers of this instance. */
};

} /* namespace partest */
#endif /* OBSERVABLE_H_ */
